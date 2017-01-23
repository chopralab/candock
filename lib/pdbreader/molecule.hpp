#ifndef MOLECULE_H
#define MOLECULE_H
#include "geom3d/geom3d.hpp"
#include "geom3d/matrix.hpp"
#include "it.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "assembly.hpp"

using namespace std;

namespace OMMIface {
	struct ForceField;
};


namespace Molib {
	class Molecules;
	class Molecule;
	class Assembly;
	class Model;
	class Chain;
	class Residue;
	class Atom;
	class Unique;
	
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
		Molecule(const Molecule &rhs) : __modified(rhs.__modified), 
			__site(rhs.__site), __name(rhs.__name), __bio_rota(rhs.__bio_rota), __bio_chain(rhs.__bio_chain) { 
			for (auto &assembly : rhs) { 
				dbgmsg("Copy constructor : molecule");
				add(new Assembly(assembly));
			}
			regenerate_bonds(rhs);
		}
		Molecule(const Molib::Molecule &rhs, const Geom3D::Point::Vec &crds);
		typedef enum {first_bio, all_bio} bio_how_many;

		string get_chain_ids(const unsigned int hm) const;


		Assembly& asym() { return this->first(); }
		
		Atom::Vec get_atoms(const string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Residue::Vec get_residues() const;

		Geom3D::Point::Vec get_crds(const string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		
		double max_dist() const;
		double max_dist(const Atom &atom) const;
		Atom& get_center_atom() const;
		void undo_mm_specific();
		void add_modified(Residue::res_tuple2 t) { __modified.insert(t); }
		bool is_modified(Residue::res_tuple2 t2) const { return (__modified.find(t2) != __modified.end()); }
		void add_site(const string &name, Residue::res_tuple2 t) { __site.insert(make_pair(name, t)); }
		void set_name(const string &name) { __name = name; }
		const string& name() const { return __name; }
		
		Molecule& filter(const unsigned int hm=Residue::notassigned, const string &chain_ids="");
		Molecule& regenerate_bonds(const Molecule&);
		Molecule& compute_overlapping_rigid_segments(Unique&);

		void prepare_for_mm(const OMMIface::ForceField &ffield, const Atom::Grid &grid);

		Geom3D::Coordinate compute_geometric_center() const { return Geom3D::compute_geometric_center(this->get_crds()); }
		double compute_rmsd(const Molecule&) const;
		double compute_rmsd_ord(const Molecule&) const;
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false);
		Assembly& add(Assembly *a) { return this->aadd(a->number(), a, this); }
		
		void init_bio(bio_how_many bhm=Molecule::all_bio);
		void add_bio_rota(const int &biomolecule_number, const int &rn, const int &row, const Geom3D::Matrix::matrix_tuple &t) {
			__bio_rota[biomolecule_number][rn].set_row(row, t); // needs copy constructor in Geom3D::Matrix
		}
		void add_bio_chain(const int &biomolecule_number, set<char> bio_chain) {
			__bio_chain[biomolecule_number] = bio_chain;
		}
		Molecule& erase_properties() { for (auto &assembly : *this) assembly.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Molecule& m);
		static string print_complex(Molecule &ligand, Molecule &receptor, const double energy, const int model = 1);
	};
	
} // Molib
#endif
	
