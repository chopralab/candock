#include "design.hpp"

#include <set>
#include <pdbreader/hydrogens.hpp>
#include <pdbreader/bondtype.hpp>

namespace design {

	Design::Design(const Molib::Molecule &start) : __original(start) {
		if (__original.empty()) {
			throw Error("Invalid molecule given for modificatiton");
		}

		// Create a new molecule, original bonds will be copied
		Molib::Hydrogens::compute_hydrogen(__original.get_atoms());
		__original.erase_properties();
	}

	const Molib::Molecules& Design::prepare_designs( const std::string& seeds_file ) {
		for ( auto &molecule : __designs ) {

			Molib::Chain   &first_chain  = molecule.first().first().first();
			Molib::Residue &main_residue = first_chain.first();
			
			auto starting_residue = first_chain.begin();
			starting_residue++;

			for ( auto residue = starting_residue; residue != first_chain.end(); ++residue )  {
				if ( residue->resi() == 1 )
					continue;
				for ( auto &atom : *residue ) {
					main_residue.add(new Molib::Atom(atom)).clear_bonds();
				}
			}

			for ( auto residue = starting_residue; residue != first_chain.end(); ++residue )  {
				if ( residue->resi() == 1 )
					continue;

				for ( auto &bond : Molib::get_bonds_in(residue->get_atoms()) ) {
					Molib::Atom &copy   = main_residue.element(bond->atom1().atom_number());
					Molib::Atom &bondee = main_residue.element(bond->atom2().atom_number());
					copy.connect( bondee ).set_bo(bond->get_bo());
				}

				//TODO: Clean this up
				for ( auto &bond : Molib::get_bonds_in(residue->get_atoms(),false) ) {
					if ( bond->atom1().br().resi() == 1 && bond->atom2().br().resi() != 1) {

						int atom_to_remove = -1;

						for ( int i = 0; i < bond->atom1().size(); ++i ) {
							if ( bond->atom2().atom_number() == bond->atom1()[i].atom_number() ) {
								atom_to_remove = i;
								break;
							}
						}

						if (atom_to_remove == -1)
							throw Error("Odd connection");

						bond->atom1().connect(main_residue.element(bond->atom2().atom_number())).set_bo(bond->get_bo());
						bond->atom1().erase(atom_to_remove);
						bond->atom1().erase_bond(bond->atom2());
						break;
					} 
					
					if ( bond->atom2().br().resi() == 1 && bond->atom1().br().resi() != 1) {

						int atom_to_remove = -1;

						for ( int i = 0; i < bond->atom2().size(); ++i ) {
							if ( bond->atom1().atom_number() == bond->atom2()[i].atom_number() ) {
								atom_to_remove = i;
								break;
							}
						}

						if (atom_to_remove == -1)
							throw Error("Odd connection");

						bond->atom2().connect(main_residue.element(bond->atom1().atom_number())).set_bo(bond->get_bo());
						bond->atom2().erase(atom_to_remove);
						bond->atom2().erase_bond(bond->atom1());
						break;
					}
				}
				
				first_chain.erase( Molib::Residue::res_pair(residue->resi(), residue->ins_code()) );
			}
		}

		return __designs.compute_hydrogen()
						.compute_bond_order()
				        .compute_bond_gaff_type()
				        .refine_idatm_type()
				        .erase_hydrogen()  // needed because refine changes connectivities
				        .compute_hydrogen()   // needed because refine changes connectivities
				        .compute_ring_type()
				        .compute_gaff_type()
				        .compute_rotatable_bonds() // relies on hydrogens being assigned
				        .erase_hydrogen()
				        .compute_overlapping_rigid_segments(seeds_file);
	}

	void Design::functionalize_hydrogens_with_fragments(const Molib::NRset& nr) {

		if (nr.size() == 0) {
			throw Error("No seeds given for modificatiton");
		}

		for ( const auto &seed : nr ) {

			// Copy seed as new residue, this will have its bonds in place
			const Molib::Residue &r = seed.first().first().first().first().first();
			Molib::Residue asmb(r);
			asmb.regenerate_bonds(r);
			asmb.set_resi(__original.first().first().first().size() + 1);
			asmb.renumber_atoms(__original.get_atoms().size());

			// "atom" is the search atom (only used for its position and type)
			for ( auto &search_atom : asmb ) {

				// Check if the atom is "exposed", I.E. not in a ring or the middle of a chain
				if (search_atom.get_bonds().size() != 1 || search_atom.element() == Molib::Element::H) {
					continue;
				}

				for ( auto &start_atom : __original.get_atoms() ) {

					// Make sure the ligand atom has an open valency
					if (start_atom->get_num_hydrogens() == 0 || start_atom->element() == Molib::Element::H) {
						continue;
					}

					// Check to see if the atom types are the same
					if (search_atom.idatm_type() == start_atom->idatm_type()) {
						Molib::Molecule modificatiton(__original);

						// Copy the new molecule into the returnable object
						__designs.add( new Molib::Molecule(modificatiton) );

						// Remove all hydrogens from the original ligand (they are not needed anymore)
						Molib::Hydrogens::erase_hydrogen(__designs.last().get_atoms());

						// Add new fragment as a "residue" of the molecule
						Molib::Chain &chain = __designs.last().first().first().first();
						chain.add(new Molib::Residue(asmb));
						chain.last().regenerate_bonds(asmb);

						// Create the relevent bond between the ligand and fragment, remove the "search atom"
						Molib::Atom &atom2  = chain.last ().element( search_atom.atom_number() );
						Molib::Atom &start2 = chain.first().element( start_atom->atom_number() );
						Molib::Atom &mod_atom = atom2.get_bonds().front()->second_atom(atom2);
						mod_atom.connect( start2 ).set_bo(1);

						__designs.last().set_name( __original.name() +"_design_with_" + seed.name() + 
							"_on_" + std::to_string(start_atom->atom_number()) +
							"_to_" + std::to_string(mod_atom.atom_number()) );

						// Fix internal atom graph
						for (int i = 0; i < mod_atom.size(); ++i) {
							if (mod_atom[i].atom_number() == atom2.atom_number()) {
								mod_atom.erase(i);
								break;
							}
						}

						mod_atom.erase_bond(atom2);
						chain.last().erase(atom2.atom_number());
					}
				}
			}
		}
	}

	void Design::functionalize_hydrogens_with_single_atoms(const std::vector<std::string>& idatms)
	{
		if ( idatms.empty() )
			throw Error("No atom types given for addition to the molecule");

		for ( auto &atom_type  : idatms ) {
			for ( auto &start_atom : __original.get_atoms() ) {
				// Make sure the ligand atom has an open valency
				if (start_atom->get_num_hydrogens() == 0 || start_atom->element() == Molib::Element::H) {
					continue;
				}

				// Copy the new molecule into the returnable object
				__designs.add( new Molib::Molecule(__original) );
				__designs.last().set_name( __original.name() + "_added_" + atom_type + "_on_" + std::to_string(start_atom->atom_number()));

				// Remove all hydrogens from the original ligand (they are not needed anymore)
				Molib::BondOrder::compute_bond_order(__designs.last().get_atoms());
				Molib::Hydrogens::erase_hydrogen(    __designs.last().get_atoms());

				Molib::Residue& mod_residue = __designs.last().first().first().first().first();
				Molib::Atom& new_atom = ( atom_type == "C" || atom_type == "O" || atom_type == "N" || atom_type == "S" ) ?
					mod_residue.add(new Molib::Atom(mod_residue.size() + 1, atom_type, Geom3D::Coordinate(0,0,0), help::idatm_mask.at(atom_type + "3") ) ) :
					mod_residue.add(new Molib::Atom(mod_residue.size() + 1, atom_type, Geom3D::Coordinate(0,0,0), help::idatm_mask.at(atom_type ) ) ) ;

				mod_residue.element(start_atom->atom_number()).connect(new_atom).set_bo(1);
			}
		}
	}

	void Design::functionalize_extremes_with_single_atoms(const std::vector< std::string >& idatms) {
		if ( idatms.empty() )
			throw Error("No atom types given for addition to the molecule");

		for ( auto &atom_type  : idatms ) {
			for ( auto &start_atom : __original.get_atoms() ) {
				// Make sure the ligand atom has an open valency
				if (start_atom->element() == Molib::Element::H || start_atom->get_bonds().size() - start_atom->get_num_hydrogens() != 1 ) {
					continue;
				}

				// Copy the new molecule into the returnable object
				__designs.add( new Molib::Molecule(__original) );
				__designs.last().set_name( __original.name() + "_changed_" + std::to_string(start_atom->atom_number()) + "_" + atom_type );

				// Remove all hydrogens from the original ligand (they are not needed anymore)
				Molib::BondOrder::compute_bond_order(__designs.last().get_atoms());
				Molib::Hydrogens::erase_hydrogen(    __designs.last().get_atoms());

				Molib::Residue& mod_residue = __designs.last().first().first().first().first();

				Molib::Atom& atom_to_change = mod_residue.element(start_atom->atom_number());
				atom_to_change.set_atom_name(atom_type);
				atom_to_change.set_element(Molib::Element(atom_type));
				atom_type == "C" || atom_type == "O" || atom_type == "N" || atom_type == "S" ? 
					atom_to_change.set_idatm_type(atom_type + "3") :
					atom_to_change.set_idatm_type(atom_type);
				atom_to_change.set_gaff_type("???");
			}
		}

	}

}
