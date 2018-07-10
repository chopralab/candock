#ifndef MOLECULES_H
#define MOLECULES_H
#include "candock/geometry/geometry.hpp"
#include "candock/molib/it.hpp"
#include "candock/molib/grid.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/residue.hpp"
#include "candock/molib/chain.hpp"
#include "candock/molib/model.hpp"
#include "candock/molib/assembly.hpp"
#include "candock/molib/molecule.hpp"

namespace candock {

namespace molib {
	class NRset;
	
	class Molecules : public template_map_container<Molecule, Molecules, NRset> {
	std::string __name; // nr-pdb name
	public:
		Molecules() : __name("") {}
		Molecules(const std::string &name) : __name(name) {}
		Molecules(const Molecules &rhs) : __name(rhs.__name) { 
			for (auto &molecule : rhs) { dbgmsg("Copy constructor : molecules");add(new Molecule(molecule)); } 
		}
		Molecule& add(Molecule *m) { return this->aadd(m, this); }
		void add(const Molecules& rhs) { for (auto &molecule : rhs) { add(new Molecule(molecule)); } }
		void set_name(const std::string &name) { __name = name; }
		const std::string& name() const { return __name; }
		void rotate(const geometry::Matrix &rota, const bool inverse=false);
		double compute_max_radius() const;
		geometry::Coordinate compute_geometric_center() const;
		
		Molecule::Vec get_molecules(const Residue::res_type &rest) const;
		geometry::Point::Vec get_crds(const std::string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Atom::Vec get_atoms(const std::string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Residue::Vec get_residues() const;

		Molecules& compute_idatm_type();
		Molecules& compute_hydrogen();
		Molecules& compute_bond_order();
		Molecules& compute_bond_gaff_type();
                Molecules& compute_chirality();
		Molecules& refine_idatm_type();
		Molecules& erase_hydrogen();
		Molecules& compute_ring_type();
		Molecules& compute_gaff_type();
		Molecules& compute_rotatable_bonds();
		
		Molecules& compute_overlapping_rigid_segments(const std::string &seeds_file="");
		Molecules& erase_properties() { for (auto &molecule : *this) molecule.erase_properties(); return *this; }

		std::set<int> get_idatm_types(std::set<int> previous=std::set<int>()) const;
		

		friend std::ostream& operator<< (std::ostream& stream, const Molecules& m);
	};


        void create_mols_from_seeds(std::set<int> &added, molib::Molecules &seeds, const molib::Molecules &mols);

} // molib

}

#endif
