#ifndef MOLECULE_H
#define MOLECULE_H
#include "candock/geometry/geometry.hpp"
#include "candock/geometry/matrix.hpp"
#include "candock/molib/it.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/residue.hpp"
#include "candock/molib/assembly.hpp"

#include <mutex>

namespace candock {

namespace OMMIface {
	struct ForceField;
};

namespace molib {
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
		typedef std::vector<Molecule*> Vec;
		
	private:
		std::set<Residue::res_tuple2> __modified;
		std::multimap<std::string, Residue::res_tuple2> __site;
	std::string __name; // usually pdb file
		typedef std::map<int, geometry::Matrix> M0;
		typedef std::map<int, M0> M1;
		M1 __bio_rota;
		std::map<int, std::set<char> > __bio_chain;
	public:
		Molecule(const std::string name) : __name(name) {}
		Molecule(const Molecule &rhs) : __modified(rhs.__modified), 
			__site(rhs.__site), __name(rhs.__name), __bio_rota(rhs.__bio_rota), __bio_chain(rhs.__bio_chain) { 
			for (auto &assembly : rhs) { 
				dbgmsg("Copy constructor : molecule");
				add(new Assembly(assembly));
			}
			regenerate_bonds(rhs);
		}

		Molecule& operator=(const Molecule &rhs) {
			__modified = rhs.__modified;
			__site = rhs.__site;
			__name = rhs.__name;
			__bio_rota = rhs.__bio_rota;
			__bio_chain = rhs.__bio_chain;
			this->clear();
			for (auto &assembly : rhs) {
				dbgmsg("Copy constructor : molecule");
				add(new Assembly(assembly));
			}
			regenerate_bonds(rhs);
			return *this;
		}

		Molecule(const molib::Molecule &rhs, const geometry::Point::Vec &crds);
		typedef enum {first_bio, all_bio} bio_how_many;

                std::string get_chain_ids(const unsigned int hm) const;

                void change_residue_name(const std::string &resn);
                void change_residue_name(std::mutex &mtx, int &ligand_cnt);

		Assembly& asym() { return this->first(); }
		std::set<int> get_idatm_types(std::set<int> previous = std::set<int>()) const;
		Atom::Vec get_atoms(const std::string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Residue::Vec get_residues() const;

		geometry::Point::Vec get_crds(const std::string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		
		double max_dist() const;
		double max_dist(const Atom &atom) const;
		Atom& get_center_atom() const;
		void undo_mm_specific();
		void add_modified(Residue::res_tuple2 t) { __modified.insert(t); }
		bool is_modified(Residue::res_tuple2 t2) const { return (__modified.find(t2) != __modified.end()); }
		void add_site(const std::string &name, Residue::res_tuple2 t) { __site.insert(make_pair(name, t)); }
		void set_name(const std::string &name) { __name = name; }
		const std::string& name() const { return __name; }
		
		Molecule& filter(const unsigned int hm=Residue::notassigned, const std::string &chain_ids="");
		Molecule& regenerate_bonds(const Molecule&);
		Molecule& compute_overlapping_rigid_segments(Unique&);

		void prepare_for_mm(const OMMIface::ForceField &ffield, const Atom::Grid &grid);

		geometry::Coordinate compute_geometric_center() const { return geometry::compute_geometric_center(this->get_crds()); }
		double compute_rmsd(const Molecule&) const;
		double compute_rmsd_ord(const Molecule&) const;
		void rotate(const geometry::Matrix &rota, const bool inverse=false);
		Assembly& add(Assembly *a) { return this->aadd(a->number(), a, this); }
		
		void init_bio(bio_how_many bhm=Molecule::all_bio);
		void add_bio_rota(const int &biomolecule_number, const int &rn, const int &row, const geometry::Matrix::matrix_tuple &t) {
			__bio_rota[biomolecule_number][rn].set_row(row, t); // needs copy constructor in geometry::Matrix
		}
		void add_bio_chain(const int &biomolecule_number, std::set<char> bio_chain) {
			__bio_chain[biomolecule_number] = bio_chain;
		}
		Molecule& erase_properties() { for (auto &assembly : *this) assembly.erase_properties(); return *this; }
		friend std::ostream& operator<< (std::ostream& stream, const Molecule& m);
	};
	
} // molib

}

#endif
