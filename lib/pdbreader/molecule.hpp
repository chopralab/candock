#ifndef MOLECULE_H
#define MOLECULE_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "geom3d/coordinate.hpp"
#include "geom3d/matrix.hpp"
#include "fragmenter/fragmenter.hpp"
#include "helper/help.hpp"
#include "graph/graph.hpp"
#include "it.hpp"
#include "element.hpp"
#include "grid.hpp"
#include "atom.hpp"

using namespace std;

namespace OMMIface {
	struct ForceField;
};

namespace Molib {
	class NRset;
	class Molecules;
	class Molecule;
	class Assembly;
	class Model;
	class Chain;
	class Residue;
	class Atom;
	class Unique;
	
	class Residue : public template_map_container<Atom, Residue, Chain, int> {
	public:
		typedef enum {notassigned=1, protein=2, nucleic=4, ion=8, water=16, hetero=32} res_type;
		typedef tuple<char, string, int, char> res_tuple2;
		typedef pair<int, char> res_pair;
	private:
		string __resn;
		int __resi;
		char __ins_code;
		Residue::res_type __rest;
		Geom3D::Coordinate __crd; // geometric center
	public:
		Residue(string resn, int resi, char ins_code, res_type(rest)) 
			: __resn(resn), __resi(resi), __ins_code(ins_code), __rest(rest) {}
		Residue(const Residue &rhs) : __resn(rhs.__resn), __resi(rhs.__resi), 
			__ins_code(rhs.__ins_code), __rest(rhs.__rest), __crd(rhs.__crd) { 
			for (auto &atom : rhs) { 
				dbgmsg("Copy constructor : residue");
				add(new Atom(atom)); 
			} 
		}
		void init_bio(Residue &residue_asym, const Geom3D::Matrix &bio_rota) {
			for (Atom &atom : residue_asym) {
				Atom &last = add(new Atom(atom));
				last.crd().rotate_inline(bio_rota); // do the actual rotation
			}
		}
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false) {
			if (inverse) for (auto &atom : *this) {	atom.crd().inverse_rotate_inline(rota); }
			else for (auto &atom : *this) {	atom.crd().rotate_inline(rota); }
		}
		Geom3D::Coordinate& crd() { return __crd; }
		void set_crd() { for (auto &atom : *this) { __crd = __crd + atom.crd(); } 
			__crd = __crd / this->size(); } // calculate geom center
		Atom& add(Atom *a) { return this->aadd(a->atom_number(), a, this); }
		string resn() const { return __resn; }
		void set_resn(const string &resn) { __resn = resn; }
		int resi() const { return __resi; }
		char ins_code() const { return __ins_code; }
		Residue::res_type rest() const { return __rest; }
		Atom::Vec get_atoms() const;
		Residue& erase_properties() { for (auto &atom : *this) atom.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Residue& r);
	};
	
	class Chain : public template_map_container<Residue, Chain, Model, Residue::res_pair> {
		char __chain_id;
		Geom3D::Coordinate __crd; // geometric center
	public:
		Chain(char chain_id) : __chain_id(chain_id) {}
		Chain(const Chain &rhs) : __chain_id(rhs.__chain_id), __crd(rhs.__crd) { 
			for (auto &residue : rhs) { 
				dbgmsg("Copy constructor : chain");
				add(new Residue(residue)); 
			} 
		}
		void init_bio(Chain &chain_asym, const Geom3D::Matrix &bio_rota) {
			for (auto &residue : chain_asym) {
				Residue &last = add(new Residue(residue)); // pair and tuple !!
				last.init_bio(residue, bio_rota);
			}
		}
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false) {
			for (auto &residue : *this) {
				residue.rotate(rota, inverse);
			}
		}
		Geom3D::Coordinate& crd() { return __crd; }
		Residue& add(Residue *r) { return this->aadd(Residue::res_pair(r->resi(), r->ins_code()), r, this); }
		Residue& residue(Residue::res_pair p) const { return this->element(p); }
		void set_crd() { 
			int sz=0; 
			for (auto &residue : *this) { 
				residue.set_crd(); 
				if (residue.rest() == Residue::protein || residue.rest() == Residue::nucleic) { 
					__crd = __crd + residue.crd(); 
					sz++; 
				}
			}
			if (sz !=0)
				__crd = __crd / sz;} // calculate geom center
		bool has_residue(Residue::res_pair p) { return this->has_element(p); }
		char chain_id() const { return __chain_id; }
		Atom::Vec get_atoms(const Residue::res_type &rest=Residue::res_type::notassigned) const;
		Chain& erase_properties() { for (auto &residue : *this) residue.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Chain& c);
	};
	
	class Model : public template_map_container<Chain, Model, Assembly, char> {
		int __number;
		map<pair<int, Residue::res_tuple2>, map<Residue::res_tuple2, Residue::res_tuple2>> __remarks;
		Fragmenter::Fragment::Vec __rigid;
	public:
		Model(int number) : __number(number) {}
		Model(const Model &rhs) : __number(rhs.__number), __remarks(rhs.__remarks) { 
			for (auto &chain : rhs) { 
				dbgmsg("Copy constructor : model");
				add(new Chain(chain)); 
			} 
		}
		Fragmenter::Fragment::Vec& get_rigid() { return __rigid; }
		const Fragmenter::Fragment::Vec& get_rigid() const { return __rigid; }
		void set_rigid(const Fragmenter::Fragment::Vec &rigid) { __rigid = rigid; }
		void init_bio(const Model &model_asym, const Geom3D::Matrix &matrix, const set<char> &chains) {
			for (const char &chain_id : chains) { 
				if (model_asym.has_chain(chain_id)) {
					Chain &last = add(new Chain(chain_id));
					last.init_bio(model_asym.chain(chain_id), matrix);
				}
			}
		}
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false) {
			for (auto &chain : *this) {
				chain.rotate(rota, inverse);
			}
		}
		Chain& add(Chain *chain) { return this->aadd(chain->chain_id(), chain, this); }
		void set_remark(const int remark_number, const Residue::res_tuple2 &ligand, 
			pair<const Residue::res_tuple2&, const Residue::res_tuple2&> rpair) { 
				this->__remarks[make_pair(remark_number, ligand)]
					.insert(make_pair(rpair.first, rpair.second)); 
		}
		Model& regenerate_bonds(const Model&);
		void set_number(const int &i) { __number = i; }
		Chain& chain(const char chain_id) const { return this->element(chain_id); }
		bool has_chain(const char chain_id) const { return this->has_element(chain_id); }
		int number() const { return __number; }
		bool remarks(const int remark_number, const Residue::res_tuple2 ligand, 
			map<Residue::res_tuple2, Residue::res_tuple2> &r) { 
			auto i = __remarks.find(make_pair(remark_number, ligand));
			if (i == __remarks.end()) return false;
			r = i->second;
			return true;
		}
		Atom::Vec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned) const;
		Model& erase_properties() { for (auto &chain : *this) chain.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Model& m);
	};
	
	class Assembly : public template_map_container<Model, Assembly, Molecule> {
		int __number;
		string __name; // ASYMMETRIC UNIT OR BIOLOGICAL ASSEMBLY
	public:
		Assembly(int number, const string name="ASYMMETRIC UNIT") : __number(number), __name(name) {}
		Assembly(const Assembly &rhs) : __number(rhs.__number), __name(rhs.__name) { 
			for (auto &model : rhs) { 
				dbgmsg("Copy constructor : assembly");
				add(new Model(model)); 
			} 
		}
		void init_bio(const Assembly &asym, map<int, Geom3D::Matrix> &matrices, const set<char> &chains) {
			for (auto &i : matrices) {
				int matrix_number = i.first;
				Geom3D::Matrix &matrix = i.second;
				Model &last = add(new Model(matrix_number));
				last.init_bio(asym.first(), matrix, chains); // asymmetric unit has just one model - the 1-st
			}
		}
		void set_number(int number) { __number = number; }
		void set_name(const string &name) { __name = name; }
		int number() const { return __number; }
		string name() const { return __name; }
		Model& add(Model *m) { return this->aadd(m->number(), m, this); }
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false) {
			for (auto &model : *this) {
				model.rotate(rota, inverse);
			}
		}
		Atom::Vec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned,
			const int model_number=-1) const;
		Assembly& erase_properties() { for (auto &model : *this) model.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Assembly& a);
	};
	
	class Molecule : public template_map_container<Assembly, Molecule, Molecules> {
	public:
		typedef vector<Molecule*> Vec;
		
	private:
		set<Residue::res_tuple2> __modified;
		multimap<string, Residue::res_tuple2> __site;
		string __name; // usually pdb file
		typedef map<int, Geom3D::Matrix> M0;
		typedef map<int, M0> M1;
		M1 __bio_rota;
		map<int, set<char> > __bio_chain;
	public:
		Molecule(const string name) : __name(name) {}
		Molecule(const Molecule &rhs) : __name(rhs.__name), __modified(rhs.__modified), 
			__site(rhs.__site), __bio_rota(rhs.__bio_rota), __bio_chain(rhs.__bio_chain) { 
			for (auto &assembly : rhs) { 
				dbgmsg("Copy constructor : molecule");
				add(new Assembly(assembly));
			}
			regenerate_bonds(rhs); 
		}
		typedef enum {first_bio, all_bio} bio_how_many;
		Assembly& asym() { return this->first(); }
		Atom::Vec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned,
			const int model_number=-1) const;
		double max_dist() const;
		double max_dist(const Atom &atom) const;
		Atom& get_center_atom() const;
		void prepare_for_mm(const OMMIface::ForceField &ff, Atom::Grid &grid);
		void undo_mm_specific();
		void add_modified(Residue::res_tuple2 t) { __modified.insert(t); }
		bool is_modified(Residue::res_tuple2 t2) const { return (__modified.find(t2) != __modified.end()); }
		void add_site(const string &name, Residue::res_tuple2 t) { __site.insert(make_pair(name, t)); }
		void set_name(const string &name) { __name = name; }
		const string& name() const { return __name; }
		Molecule& filter(const unsigned int hm=Residue::notassigned, const string &chain_ids="");
		Molecule& compute_idatm_type();
		Molecule& refine_idatm_type();
		Molecule& compute_gaff_type();
		Molecule& compute_bond_order();
		Molecule& compute_ring_type();
		Molecule& compute_bond_gaff_type();
		Molecule& compute_hydrogen();
		Molecule& erase_hydrogen();
		Molecule& regenerate_bonds(const Molecule&);
		Molecule& compute_rotatable_bonds();
		Molecule& compute_overlapping_rigid_segments(Unique&);
		Molecule& compute_seeds(Unique &u);

		Geom3D::Coordinate compute_geometric_center() const { 
			Geom3D::Coordinate center;
			for (auto &patom : this->get_atoms()) {
				center = center + patom->crd();
			}
			center = center / this->get_atoms().size();
			return center;
		}
		double compute_rmsd(const Molecule&) const;
		double compute_rmsd_ord(const Molecule&) const;
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false) {
			for (auto &assembly : *this) {
				assembly.rotate(rota, inverse);
			}
		}
		Assembly& add(Assembly *a) { return this->aadd(a->number(), a, this); }
		void init_bio(bio_how_many bhm=Molecule::all_bio) {
			for(auto &i : __bio_rota) {
				int biomolecule_number = i.first;
				M0 &matrices = i.second;
				Assembly &last = add(new Assembly(biomolecule_number, "BIOLOGICAL ASSEMBLY"));
				last.init_bio(asym(), matrices, __bio_chain[biomolecule_number]);
				if (bhm == Molecule::first_bio) { // exit after initializing the first assembly
					break;
				}
			}
		}
		void add_bio_rota(const int &biomolecule_number, const int &rn, const int &row, const Geom3D::Matrix::matrix_tuple &t) {
			__bio_rota[biomolecule_number][rn].set_row(row, t); // needs copy constructor in Geom3D::Matrix
		}
		void add_bio_chain(const int &biomolecule_number, set<char> bio_chain) {
			__bio_chain[biomolecule_number] = bio_chain;
		}
		Molecule& erase_properties() { for (auto &assembly : *this) assembly.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Molecule& m);
	};
	
	class Molecules : public template_map_container<Molecule, Molecules, NRset> {
		string __name; // nr-pdb name
	public:
		Molecules() : __name("") {}
		Molecules(const string &name) : __name(name) {}
		Molecules(const Molecules &rhs) : __name(rhs.__name) { 
			for (auto &molecule : rhs) { dbgmsg("Copy constructor : molecules");add(new Molecule(molecule)); } 
		}
		Molecule& add(Molecule *m) { return this->aadd(m, this); }
		void set_name(const string &name) { __name = name; }
		const string& name() const { return __name; }
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false) {
			for (auto &molecule : *this) {
				molecule.rotate(rota, inverse);
			}
		}
		double compute_max_radius() const { 
			Geom3D::Coordinate gc = compute_geometric_center();
			double max_radius = -HUGE_VAL;
			for (auto &molecule : *this) {
				Atom::Vec atoms = molecule.get_atoms();
				for (auto &patom : atoms)
					max_radius = max(max_radius, patom->crd().distance(gc));
			}
			return max_radius;
		}
		Geom3D::Coordinate compute_geometric_center() const { 
			Geom3D::Coordinate center;
			for (auto &molecule : *this) {
				center = center + molecule.compute_geometric_center();
			}
			return center / this->size();
		}
		Atom::Vec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned,
			const int model_number=-1) const;
		Molecules& compute_ring_type() { for (auto &molecule : *this) 
			molecule.compute_ring_type(); return *this; }
		Molecules& compute_bond_gaff_type() { for (auto &molecule : *this) 
			molecule.compute_bond_gaff_type(); return *this; }
		Molecules& compute_idatm_type() { for (auto &molecule : *this) 
			molecule.compute_idatm_type(); return *this; }
		Molecules& refine_idatm_type() { for (auto &molecule : *this) 
			molecule.refine_idatm_type(); return *this; }
		Molecules& compute_gaff_type() { for (auto &molecule : *this) 
			molecule.compute_gaff_type(); return *this; }
		Molecules& compute_bond_order() { for (auto &molecule : *this) 
			molecule.compute_bond_order(); return *this; }
		Molecules& compute_hydrogen() { for (auto &molecule : *this) 
			molecule.compute_hydrogen(); return *this; }
		Molecules& erase_hydrogen() { for (auto &molecule : *this) 
			molecule.erase_hydrogen(); return *this; }
		Molecules& compute_rotatable_bonds() { for (auto &molecule : *this) 
			molecule.compute_rotatable_bonds(); return *this; }
		Molecules& compute_overlapping_rigid_segments(const string &seeds_file="");
		Molecules& erase_properties() { for (auto &molecule : *this) molecule.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Molecules& m);
	};
	
	class NRset : public template_map_container<Molecules, NRset, NRset> {
	public:
		Molecules& add(Molecules *m) { return this->aadd(m, this); }
		NRset& erase_properties() { for (auto &molecules : *this) molecules.erase_properties(); return *this; }
		Atom::Vec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned,
			const int model_number=-1) const;
		friend ostream& operator<< (ostream& stream, const NRset& m);
	};
	
	
	set<int> get_idatm_types(const Molecules&, set<int> previous=set<int>());


} // Molib
#endif
	
