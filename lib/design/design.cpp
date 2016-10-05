#include "design.hpp"

#include <set>
#include <pdbreader/hydrogens.hpp>
#include <pdbreader/bondtype.hpp>

namespace design {
	Molib::Molecules Design::add_fragments_to_existing_molecule(const Molib::Molecule& start, const Molib::NRset& nr) {

		const string original_name = start.name();

		if (nr.size() == 0) {
			throw Error("No seeds given for modificatiton");
		}
		
		if (start.empty()) {
			throw Error("Invalid molecule given for modificatiton");
		}

		Molib::Molecules all_designs;

		for ( const auto &seed : nr ) {

			cout << "Attempting seed " << seed.name() << endl;

			// Create a new molecule, original bonds will be copied
			Molib::Molecule modificatiton(start);
			Molib::Hydrogens::compute_hydrogen(modificatiton.get_atoms());
			Molib::BondOrder::compute_bond_order(modificatiton.get_atoms());

			// Copy seed as new residue, this will have its bonds in place
			const Molib::Residue &r = seed.first().first().first().first().first();
			Molib::Residue asmb(r);
			asmb.regenerate_bonds(r);
			asmb.set_resi(modificatiton.first().first().first().size() + 1);
			asmb.renumber_atoms(modificatiton.get_atoms().size());

			// "atom" is the search atom (only used for its position and type)
			for ( auto &atom : asmb ) {

				// Check if the atom is "exposed", I.E. not in a ring or the middle of a chain
				if (atom.get_bonds().size() != 1 || atom.element() == Molib::Element::H) {
					continue;
				}

				for ( auto &start_atom : modificatiton.get_atoms() ) {

					// Make sure the ligand atom has an open valency
					if (start_atom->get_num_hydrogens() == 0 || start_atom->element() == Molib::Element::H) {
						continue;
					}

					// Check to see if the atom types are the same
					if (atom.idatm_type_unmask() == start_atom->idatm_type_unmask()) {
						modificatiton.set_name( original_name +"_design_with_" + seed.name() + "_on_" + std::to_string(start_atom->atom_number()) );

						// Copy the new molecule into the returnable object
						all_designs.add( new Molib::Molecule(modificatiton) );

						// Remove all hydrogens from the original ligand (they are not needed anymore)
						Molib::Hydrogens::erase_hydrogen(all_designs.last().get_atoms());

						// Add new fragment as a "residue" of the molecule
						Molib::Chain &chain = all_designs.last().first().first().first();
						chain.add(new Molib::Residue(asmb));
						chain.last().regenerate_bonds(asmb);

						// Create the relevent bond between the ligand and fragment, remove the "search atom"
						Molib::Atom &atom2  = chain.last ().element( atom.atom_number() );
						Molib::Atom &start2 = chain.first().element( start_atom->atom_number() );
						Molib::Atom &mod_atom = atom2.get_bonds().front()->second_atom(atom2);
						mod_atom.connect( start2 );
						
						chain.last().erase(atom2.atom_number());
						mod_atom.erase_bond(atom2);
						
						Molib::Atom::Vec remover = chain.last().get_atoms();
						remover.erase( std::remove_if( remover.begin(), remover.end(), 
							[&atom2] (Molib::Atom *patom) -> bool {
								return atom2.atom_number() == patom->atom_number();
							} ), remover.end());
					}
				}
			}
		}

		return all_designs.compute_hydrogen()
				   .compute_bond_order()
				   .compute_bond_gaff_type()
				   .refine_idatm_type()
				   .erase_hydrogen()  // needed because refine changes connectivities
				   .compute_hydrogen()   // needed because refine changes connectivities
				   .compute_ring_type()
				   .compute_gaff_type()
				   .compute_rotatable_bonds() // relies on hydrogens being assigned
				   .erase_hydrogen();;
	}

}
