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
#include "it.hpp"
#include "element.hpp"
#include "fragmenter/fragmenter.hpp"
#include "helper/help.hpp"
//~ #include "fragmenter/unique.hpp"
#include "grid.hpp"
#include "bond.hpp"
#include "graph/graph.hpp"

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
	
	typedef Grid<Atom> MolGrid;
	typedef set<Atom*> AtomSet;
	typedef vector<Atom*> AtomVec;
	typedef pair<Atom*, Atom*> AtomPair;

	//~ class Atom : public template_map_container<Bond, Atom, Residue, Atom*> {
	class Atom : public template_vector_container<Atom*, Atom> {
		int __atom_number;
		string __atom_name;
		Geom3D::Coordinate __crd;
		int __idatm_type;
		string __gaff_type;
		Element __element;
		string __smiles_label;
		map<string, int> __smiles_prop;
		map<int, int> __aps;
		void *__br; // back reference
		map<Atom*, shared_ptr<Bond>> __bonds;
	public:
		Atom(const Atom &rhs) : __atom_number(rhs.__atom_number), 
			__atom_name(rhs.__atom_name), __crd(rhs.__crd), 
			__idatm_type(rhs.__idatm_type), __gaff_type(rhs.__gaff_type),
			__element(rhs.__element), __smiles_label(rhs.__smiles_label),
			__smiles_prop(rhs.__smiles_prop), __aps(rhs.__aps) {
			dbgmsg("Copy constructor : atom");
		}
		typedef tuple<int, string, unique_ptr<Geom3D::Coordinate>, double> atom_tuple;
		//~ Atom(const int atom_number, const string &smiles_label, const set<int> num_bonds, 
			//~ const set<string> props) : __atom_number(atom_number), 
			//~ __smiles_label(smiles_label), __num_bonds(num_bonds), __props(props) {}
		Atom(const int atom_number, const string &smiles_label, 
			const map<string, int> smiles_prop) : __atom_number(atom_number), 
			__smiles_label(smiles_label), __smiles_prop(smiles_prop), __element(""),
			__br(nullptr) {}
		Atom(const Geom3D::Coordinate &crd) : __crd(crd), __element("") {}
		Atom(const Geom3D::Coordinate &crd, const int &idatm_type) : __crd(crd), 
			__idatm_type(idatm_type), __element("") {}
		Atom(int atom_number, string atom_name, const Geom3D::Coordinate &crd, 
			const int idatm_type) : __atom_number(atom_number), 
			__atom_name(atom_name), __crd(crd), __idatm_type(idatm_type), 
			__gaff_type("???"), __element(atom_name) {}
			
		//~ Bond& add(Bond *b) { return this->aadd(&b->second_atom(), b, this); }
		//~ Bond& add(Bond *b) { return this->aadd_no_br(&b->second_atom(), b); }
		//~ Bond& get_bond(const Atom &atom2) const { return this->template_map_container<Bond, Atom, Residue, Atom*>::element(&atom2); }
		//~ Bond& get_bond(Atom &atom2) const { return this->
			//~ template_map_container<Bond, Atom, Residue, Atom*>::element(&atom2); }
		void clear_bonds() { __bonds.clear(); }
		BondVec get_bonds() const { BondVec bonds; for (auto &kv : __bonds) 
			bonds.push_back(&*kv.second); return bonds; }
		const shared_ptr<Bond>& get_shared_ptr_bond(Atom &other) const { return __bonds.at(&other); }
		Bond& get_bond(Atom &other) const { return *get_shared_ptr_bond(other); }
		shared_ptr<Bond>& insert_bond(Atom &other, const shared_ptr<Bond> &bond) { return __bonds.insert({&other, bond}).first->second; }
		shared_ptr<Bond>& insert_bond(Atom &other, Bond *bond) { return __bonds.insert({&other, shared_ptr<Bond>(bond)}).first->second; }
		void erase_bond(Atom &other) { __bonds.erase(&other); }
		//~ void insert_bond(Atom &other, shared_ptr<Bond> &bond) { bonds.insert({&other, bond}); }
		//~ bool has_bond(Atom &other) const { return __bonds.count(&other); }
		bool is_adjacent(Atom &other) const { return __bonds.count(&other); }
		bool is_adjacent(const string &atom_name) const { for (auto &other : *this) if (other.atom_name() == atom_name) return true; return false; }
		//~ int get_num_hydrogens() const { return __num_hydrogens; }
		//~ void set_num_hydrogens(const int num_h) { __num_hydrogens = num_h; }
		int get_num_hydrogens() const;
		//~ bool is_adjacent(Atom &atom2) { return this->has_element(&atom2); }
		int atom_number() const { return __atom_number; }
		void set_atom_name(const string &atom_name) { __atom_name = atom_name; }
		void set_idatm_type(const string &idatm_type) { __idatm_type = help::idatm_mask.at(idatm_type); }
		void set_gaff_type(const string &gaff_type) { __gaff_type = gaff_type; }
		//~ Atom& insert_property(const string &prop) { __props.insert(prop); return *this; }
		//~ Atom& erase_property(const string &prop) { __props.erase(prop); return *this; }
		//~ bool has_property(const string &prop) const { return __props.count(prop); }
		Atom& insert_property(const string &prop, const int count) { __smiles_prop.insert({prop, count}); return *this; }
		Atom& add_property(const string &prop) { __smiles_prop[prop]++; return *this; }
		Atom& erase_property(const string &prop) { __smiles_prop.erase(prop); return *this; }
		bool has_property(const string &prop) const { return __smiles_prop.count(prop); }
		int get_num_property(const string &prop) const { 
			dbgmsg("__smiles prop(" << prop << ") count = " << __smiles_prop.count(prop)); 
			if (__smiles_prop.count(prop) == 0) return 0;
			return __smiles_prop.at(prop); }
		int get_num_bond_with_bond_gaff_type(const string &prop) const;
		int compute_num_property(const string &prop) const;
		void set_crd(const Geom3D::Coordinate &crd) { __crd = crd; }
		string idatm_type_unmask() const { return help::idatm_unmask[__idatm_type]; }
		int idatm_type() const { return __idatm_type; }
		const string& gaff_type() const { return __gaff_type; }
		const string& atom_name() const { return __atom_name; }
		Element element() const { return __element; }
		Geom3D::Coordinate& crd() { return __crd; }
		const Geom3D::Coordinate& crd() const { return __crd; }
		friend ostream& operator<< (ostream& stream, const Atom& a);
		double distance() const { return 0.0; } // just dummy : needed by grid
		void distance(double d) const {} // just dummy : needed by grid
		const map<int, int>& get_aps() const { return __aps; }
		void set_members(const string &str);
		const Residue& br() const { return *static_cast<const Residue*>(__br); }
		void set_br(void *br) { __br = br; }
		// the following are required for MolGraph :-)
		bool compatible(const Atom &atom) const;
		//~ string get_label() const { return help::idatm_unmask[__idatm_type]; }
		string get_label() const { return (__smiles_label.empty() 
			? help::idatm_unmask[__idatm_type] : __smiles_label); }
		//~ string print() const { return get_label(); }
		const int weight() const { return 0; } // dummy for graph ostream operator
	};
	
	class Residue : public template_map_container<Atom, Residue, Chain, int> {
	public:
		typedef enum {protein, nucleic, ion, water, hetero, notassigned} res_type;
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
		AtomVec get_atoms() const;
		//~ BondVec get_bonds() const { BondVec bonds; for (auto &atom : *this) 
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
		AtomVec get_atoms(const Residue::res_type &rest=Residue::res_type::notassigned) const;
		friend ostream& operator<< (ostream& stream, const Chain& c);
	};
	
	class Model : public template_map_container<Chain, Model, Assembly, char> {
		int __number;
		map<pair<int, Residue::res_tuple2>, map<Residue::res_tuple2, Residue::res_tuple2>> __remarks;
		Fragments __seeds, __rigid;
	public:
		Model(int number) : __number(number) {}
		Model(const Model &rhs) : __number(rhs.__number), __remarks(rhs.__remarks) { 
			for (auto &chain : rhs) { 
				dbgmsg("Copy constructor : model");
				add(new Chain(chain)); 
			} 
		}
		Fragments& get_seeds() { return __seeds; }
		Fragments& get_rigid() { return __rigid; }
		const Fragments& get_seeds() const { return __seeds; }
		const Fragments& get_rigid() const { return __rigid; }
		void set_seeds(const Fragments &seeds) { __seeds = seeds; }
		void set_rigid(const Fragments &rigid) { __rigid = rigid; }
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
		AtomVec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned) const;
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
				add(new Model(matrix_number));
				last().init_bio(asym[0], matrix, chains); // asymmetric unit has just one model - the 0-th
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
		AtomVec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned,
			const int model_number=-1) const;
		friend ostream& operator<< (ostream& stream, const Assembly& a);
	};
	
	class Molecule : public template_map_container<Assembly, Molecule, Molecules> {
		set<Residue::res_tuple2> __modified;
		multimap<string, Residue::res_tuple2> __site;
		string __name; // usually pdb file
		typedef map<int, Geom3D::Matrix> M0;
		typedef map<int, M0> M1;
		M1 __bio_rota;
		map<int, set<char> > __bio_chain;
		//~ Fragments __seeds, __rigid;
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
		AtomVec get_atoms(const string &chain_ids="", 
			const Residue::res_type &rest=Residue::res_type::notassigned,
			const int model_number=-1) const;
		//~ Fragments& get_seeds() { return __seeds; }
		//~ Fragments& get_rigid() { return __rigid; }
		double max_dist() const;
		void prepare_for_mm(const OMMIface::ForceField &ff, MolGrid &grid);
		//~ void prepare_for_mm(const OMMIface::ForceField&);
		void undo_mm_specific();
		void add_modified(Residue::res_tuple2 t) { __modified.insert(t); }
		bool is_modified(Residue::res_tuple2 t2) const { return (__modified.find(t2) != __modified.end()); }
		void add_site(const string &name, Residue::res_tuple2 t) { __site.insert(make_pair(name, t)); }
		void set_name(const string &name) { __name = name; }
		const string& name() const { return __name; }
		//~ bool has_hydrogen() const;
		Molecule& compute_idatm_type();
		Molecule& compute_gaff_type();
		Molecule& compute_bond_order();
		Molecule& compute_ring_type();
		Molecule& compute_bond_gaff_type();
		Molecule& compute_hydrogen();
		Molecule& erase_hydrogen();
		//~ Molecule& compute_fragment(Unique&);
		Molecule& regenerate_bonds(const Molecule&);
		Molecule& compute_rotatable_bonds();
		Molecule& compute_overlapping_rigid_segments();
		Molecule& compute_seeds(Unique &u);
		double compute_rmsd(const Molecule&) const;
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
				add(new Assembly(biomolecule_number, "BIOLOGICAL ASSEMBLY"));
				last().init_bio(asym(), matrices, __bio_chain[biomolecule_number]);
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
				AtomVec atoms = molecule.get_atoms();
				for (auto &patom : atoms)
					max_radius = max(max_radius, patom->crd().distance(gc));
			}
			return max_radius;
		}
		Geom3D::Coordinate compute_geometric_center() const { 
			Geom3D::Coordinate gc;
			int i = 0;
			for (auto &molecule : *this) {
				AtomVec atoms = molecule.get_atoms();
				for (auto &patom : atoms)
					gc = gc + patom->crd();
				i += atoms.size();
			}
			return gc / (double) i;
		}
		Molecules& compute_ring_type() { for (auto &molecule : *this) 
			molecule.compute_ring_type(); return *this; }
		Molecules& compute_bond_gaff_type() { for (auto &molecule : *this) 
			molecule.compute_bond_gaff_type(); return *this; }
		Molecules& compute_idatm_type() { for (auto &molecule : *this) 
			molecule.compute_idatm_type(); return *this; }
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
		Molecules& compute_overlapping_rigid_segments() { for (auto &molecule : *this) 
			molecule.compute_overlapping_rigid_segments(); return *this; }
		//~ Molecules& compute_fragment(const string &seeds_file="") { 
			//~ Unique u(seeds_file); for (auto &molecule : *this) 
			//~ molecule.compute_fragment(u); return *this; }
		Molecules& compute_seeds(const string &seeds_file="");
			//~ Unique u(seeds_file); for (auto &molecule : *this) 
			//~ molecule.compute_seeds(u); return *this; }
		friend ostream& operator<< (ostream& stream, const Molecules& m);
	};
	
	class NRset : public template_map_container<Molecules, NRset, NRset> {
	public:
		Molecules& add(Molecules *m) { return this->aadd(m, this); }
		friend ostream& operator<< (ostream& stream, const NRset& m);
	};
	
	typedef Glib::Graph<Atom> MolGraph;

	
	double compute_rmsd(const MolGraph&, const MolGraph&, const MolGraph::Matches&);
	set<int> get_idatm_types(const Molecules&, set<int> previous=set<int>());

	//~ template<typename T>
	//~ BondVec get_bonds_in(const T &atoms, bool in=true) {
		//~ set<AtomPair> visited;
		//~ BondVec bonds;
		//~ for (auto &patom : atoms) {
			//~ Atom &atom1 = *patom;
			//~ for (auto &bond : atom1) {
				//~ Atom &atom2 = bond.second_atom();
				//~ if ((in || !atoms.count(&atom2)) 
					//~ && (!in || atoms.count(&atom2)) 
					//~ && !visited.count({&atom2, &atom1})) {
					//~ visited.insert({&atom1, &atom2});
					//~ bonds.push_back(&bond);
				//~ }
			//~ }
		//~ }
		//~ return bonds;
	//~ }
	//~ BondVec get_bonds_in(const AtomSet &atoms, bool in=true);
	//~ BondVec get_bonds_in(const AtomVec &atoms, bool in=true);
	BondSet get_bonds_in(const AtomSet &atoms, bool in=true);
	BondSet get_bonds_in(const AtomVec &atoms, bool in=true);

	MolGraph create_graph(const AtomVec &atoms);
	MolGraph create_graph(const AtomSet &atoms);
	//~ MolGraph create_mol_graph(const help::smiles &edges);

	ostream& operator<< (ostream& stream, const AtomSet&atoms);
	ostream& operator<< (ostream& stream, const AtomVec&atoms);

} // Molib
#endif
	
