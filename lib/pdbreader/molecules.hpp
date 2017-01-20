#ifndef MOLECULES_H
#define MOLECULES_H
#include "geom3d/geom3d.hpp"
#include "helper/help.hpp"
#include "it.hpp"
#include "grid.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
#include "model.hpp"
#include "assembly.hpp"
#include "molecule.hpp"

using namespace std;

namespace Molib {
	class NRset;
	
	class Molecules : public template_map_container<Molecule, Molecules, NRset> {
		string __name; // nr-pdb name
	public:
		Molecules() : __name("") {}
		Molecules(const string &name) : __name(name) {}
		Molecules(const Molecules &rhs) : __name(rhs.__name) { 
			for (auto &molecule : rhs) { dbgmsg("Copy constructor : molecules");add(new Molecule(molecule)); } 
		}
		Molecule& add(Molecule *m) { return this->aadd(m, this); }
		void add(const Molecules& rhs) { for (auto &molecule : rhs) { add(new Molecule(molecule)); } }
		void set_name(const string &name) { __name = name; }
		const string& name() const { return __name; }
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false);
		double compute_max_radius() const;
		Geom3D::Coordinate compute_geometric_center() const;
		
		Molecule::Vec get_molecules(const Residue::res_type &rest) const;
		Geom3D::Point::Vec get_crds(const string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Atom::Vec get_atoms(const string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Residue::Vec get_residues() const;

		Molecules& compute_idatm_type();
		Molecules& compute_hydrogen();
		Molecules& compute_bond_order();
		Molecules& compute_bond_gaff_type();
		Molecules& refine_idatm_type();
		Molecules& erase_hydrogen();
		Molecules& compute_ring_type();
		Molecules& compute_gaff_type();
		Molecules& compute_rotatable_bonds();
		
		Molecules& compute_overlapping_rigid_segments(const string &seeds_file="");
		Molecules& erase_properties() { for (auto &molecule : *this) molecule.erase_properties(); return *this; }

		set<int> get_idatm_types(set<int> previous=set<int>()) const;
		

		friend ostream& operator<< (ostream& stream, const Molecules& m);
	};


} // Molib
#endif
	
