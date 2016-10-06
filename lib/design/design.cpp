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
	}

	const Molib::Molecules& Design::get_prepared_designs() {
		cout << "ASDF" << endl;
		for ( auto &molecule : __designs ) {

			Molib::Chain   &first_chain  = molecule.first().first().first();
			Molib::Residue &main_residue = first_chain.first();

			auto starting_residue = first_chain.begin();
			starting_residue++;

			for ( auto residue = starting_residue; residue != first_chain.end(); ++residue )  {
				if ( residue->resi() == 1 )
					continue;
				for ( auto &atom : *residue ) {
					main_residue.add(new Molib::Atom(atom));
				}
			}
			
			for ( auto residue = starting_residue; residue != first_chain.end(); ++residue )  {
				if ( residue->resi() == 1 )
					continue;
				for ( auto &atom : *residue ) {
					for ( auto &bond : atom.get_bonds()) {
						Molib::Atom &copy = main_residue.element(atom.atom_number());
						Molib::Atom &original_bondee = bond->second_atom(atom);
						copy.connect( main_residue.element(original_bondee.atom_number()));
					}
				}
			}
			
			for ( auto residue = starting_residue; residue != first_chain.end(); ++residue )  {
				if ( residue->resi() == 1 )
					continue;
				cout << residue->resi() << " " << residue->ins_code() << endl;
				first_chain.erase( Molib::Residue::res_pair(residue->resi(), residue->ins_code()) );
			}
		}
		
		return __designs;
	}

	void Design::add_fragments_to_existing_molecule(const Molib::NRset& nr) {

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
			for ( auto &atom : asmb ) {

				// Check if the atom is "exposed", I.E. not in a ring or the middle of a chain
				if (atom.get_bonds().size() != 1 || atom.element() == Molib::Element::H) {
					continue;
				}

				for ( auto &start_atom : __original.get_atoms() ) {

					// Make sure the ligand atom has an open valency
					if (start_atom->get_num_hydrogens() == 0 || start_atom->element() == Molib::Element::H) {
						continue;
					}

					// Check to see if the atom types are the same
					if (atom.idatm_type_unmask() == start_atom->idatm_type_unmask()) {

						Molib::Molecule modificatiton(__original);

						modificatiton.set_name( __original.name() +"_design_with_" + seed.name() + "_on_" + std::to_string(start_atom->atom_number()) );

						// Copy the new molecule into the returnable object
						__designs.add( new Molib::Molecule(modificatiton) );

						// Remove all hydrogens from the original ligand (they are not needed anymore)
						Molib::Hydrogens::erase_hydrogen(__designs.last().get_atoms());

						// Add new fragment as a "residue" of the molecule
						Molib::Chain &chain = __designs.last().first().first().first();
						chain.add(new Molib::Residue(asmb));
						chain.last().regenerate_bonds(asmb);

						// Create the relevent bond between the ligand and fragment, remove the "search atom"
						Molib::Atom &atom2  = chain.last ().element( atom.atom_number() );
						Molib::Atom &start2 = chain.first().element( start_atom->atom_number() );
						Molib::Atom &mod_atom = atom2.get_bonds().front()->second_atom(atom2);
						mod_atom.connect( start2 );

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

}
