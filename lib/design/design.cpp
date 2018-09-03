#include "candock/design/design.hpp"

#include <set>

#include "statchem/fragmenter/unique.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/molib/atomtype.hpp"
#include "statchem/molib/bondtype.hpp"

using namespace std;
using namespace statchem;

namespace candock {

namespace design {

Design::Design(const molib::Molecule& start, molib::Unique& existing)
    : __original(start), __existing(existing) {
    if (__original.empty()) {
        throw Error("Invalid molecule given for modificatiton");
    }

    // Create a new molecule, original bonds will be copied
    __original.first().first().first().first().compute_hydrogen();
    __original.erase_properties();
    __original.first().first().first().first().renumber_atoms(1);
}

const molib::Molecules& Design::prepare_designs() {
    for (auto& molecule : __designs) {
        molib::Chain& first_chain = molecule.first().first().first();
        molib::Residue& main_residue = first_chain.first();

        auto starting_residue = first_chain.begin();
        starting_residue++;

        for (auto residue = starting_residue; residue != first_chain.end();
             ++residue) {
            if (residue->resi() == 1) continue;
            for (auto& atom : *residue) {
                main_residue.add(new molib::Atom(atom)).clear_bonds();
            }
        }

        for (auto residue = starting_residue; residue != first_chain.end();
             ++residue) {
            if (residue->resi() == 1) continue;

            for (auto& bond : molib::get_bonds_in(residue->get_atoms())) {
                molib::Atom& copy =
                    main_residue.element(bond->atom1().atom_number());
                molib::Atom& bondee =
                    main_residue.element(bond->atom2().atom_number());
                copy.connect(bondee).set_bo(bond->get_bo());
            }

            // TODO: Clean this up
            for (auto& bond :
                 molib::get_bonds_in(residue->get_atoms(), false)) {
                if (bond->atom1().br().resi() == 1 &&
                    bond->atom2().br().resi() != 1) {
                    int atom_to_remove = -1;

                    for (size_t i = 0; i < bond->atom1().size(); ++i) {
                        if (bond->atom2().atom_number() ==
                            bond->atom1()[i].atom_number()) {
                            atom_to_remove = i;
                            break;
                        }
                    }

                    if (atom_to_remove == -1) throw Error("Odd connection");

                    bond->atom1()
                        .connect(
                            main_residue.element(bond->atom2().atom_number()))
                        .set_bo(bond->get_bo());
                    bond->atom1().erase(atom_to_remove);
                    bond->atom1().erase_bond(bond->atom2());
                    break;
                }

                if (bond->atom2().br().resi() == 1 &&
                    bond->atom1().br().resi() != 1) {
                    int atom_to_remove = -1;

                    for (size_t i = 0; i < bond->atom2().size(); ++i) {
                        if (bond->atom1().atom_number() ==
                            bond->atom2()[i].atom_number()) {
                            atom_to_remove = i;
                            break;
                        }
                    }

                    if (atom_to_remove == -1) throw Error("Odd connection");

                    bond->atom2()
                        .connect(
                            main_residue.element(bond->atom1().atom_number()))
                        .set_bo(bond->get_bo());
                    bond->atom2().erase(atom_to_remove);
                    bond->atom2().erase_bond(bond->atom1());
                    break;
                }
            }

            first_chain.erase(
                molib::Residue::res_pair(residue->resi(), residue->ins_code()));
        }

        // Find terminal atoms that have invalid atom types for terminal atoms.
        // These checks need to be changed (improved) when dealing with correct
        // bond typing
        for (auto atom : molecule.get_atoms()) {
            if ((atom->idatm_type_unmask() == "Car" ||
                 atom->idatm_type_unmask() == "C2") &&
                atom->get_bonds().size() == 1 &&
                atom->get_bonds()
                        .at(0)
                        ->second_atom(*atom)
                        .idatm_type_unmask()
                        .back() != '2') {
                atom->set_idatm_type("C3");
            } else if ((atom->idatm_type_unmask() == "Nar" ||
                        atom->idatm_type_unmask() == "N2") &&
                       atom->get_bonds().size() == 1 &&
                       atom->get_bonds()
                               .at(0)
                               ->second_atom(*atom)
                               .idatm_type_unmask()
                               .back() != '2') {
                atom->set_idatm_type("N3");
            }
        }
    }

    return __designs;
}

bool check_clash_for_design(molib::Atom::Vec molecule, molib::Atom::Vec seed,
                            double clash_coeff, int skip) {
    for (const auto i : molecule) {
        const double vdw1 = i->radius();
        for (const auto j : seed) {
            if (j->atom_number() == skip) {
                continue;
            }
            const double vdw2 = j->radius();
            if (i->crd().distance_sq(j->crd()) <
                pow(clash_coeff * (vdw1 + vdw2), 2))
                return true;
        }
    }
    return false;
}

void design::Design::functionalize_hydrogens_with_fragments(
    const molib::NRset& nr, const double cutoff, const double clash_coeff,
    const std::tuple<double, size_t, size_t, size_t>& lipinski_values) {
    if (nr.size() == 0) {
        throw Error("No seeds given for modificatiton");
    }

    std::tuple<double, size_t, size_t, size_t> original_lipinski =
        molib::AtomType::determine_lipinski(__original.get_atoms());

    for (auto& start_atom : __original.get_atoms()) {
        // Make sure the ligand atom has an open valency
        if (start_atom->get_num_hydrogens() == 0 ||
            start_atom->element() == molib::Element::H) {
            continue;
        }

        for (const auto& fragment : nr) {
            fragment.first().first().first().first().first().compute_hydrogen();
            std::tuple<double, size_t, size_t, size_t> frag_lipinski =
                molib::AtomType::determine_lipinski(fragment.first()
                                                        .first()
                                                        .first()
                                                        .first()
                                                        .first()
                                                        .get_atoms());
            fragment.first().first().first().first().first().erase_hydrogen();

            if (std::get<0>(original_lipinski) + std::get<0>(frag_lipinski) >
                    std::get<0>(lipinski_values) ||
                std::get<1>(original_lipinski) + std::get<1>(frag_lipinski) >
                    std::get<1>(lipinski_values) ||
                std::get<2>(original_lipinski) + std::get<2>(frag_lipinski) >
                    std::get<2>(lipinski_values) ||
                std::get<3>(original_lipinski) + std::get<3>(frag_lipinski) >
                    std::get<3>(lipinski_values)) {
                continue;
            }

            std::set<std::pair<int, int> > already_added;

            // This variable controls if we've already added the seed
            for (const auto& seed : fragment) {
                // Copy seed as new residue, this will have its bonds in place
                const molib::Residue& r = seed.first().first().first().first();

                // "atom" is the search atom (only used for its position and
                // type)
                for (auto& search_atom : r) {
                    std::pair<int, int> test_already_added(
                        search_atom.idatm_type(), start_atom->idatm_type());

                    if (already_added.count(test_already_added) != 0) {
                        continue;
                    }

                    // Check if the atom is "exposed", I.E. not in a ring or the
                    // middle of a chain
                    if (search_atom.get_bonds().size() != 1 ||
                        search_atom.element() == molib::Element::H) {
                        continue;
                    }

                    // Check to see if the atom types are the same
                    if (search_atom.idatm_type() != start_atom->idatm_type()) {
                        continue;
                    }

                    // Ensure the fragment is close enough to the original
                    // ligand in the pocket
                    if (start_atom->crd().distance(search_atom.crd()) >
                        cutoff) {
                        continue;
                    }

                    if (check_clash_for_design(__original.get_atoms(),
                                               seed.get_atoms(), clash_coeff,
                                               search_atom.atom_number())) {
                        continue;
                    }

                    already_added.insert(test_already_added);

                    // Copy the new molecule into the returnable object
                    molib::Molecule& new_design =
                        __designs.add(new molib::Molecule(__original));

                    // Remove all hydrogens from the original ligand (they are
                    // not needed anymore)
                    new_design.first().first().first().first().erase_hydrogen();

                    // Add new fragment as a "residue" of the molecule
                    molib::Chain& chain = new_design.first().first().first();
                    molib::Residue* new_res = new molib::Residue(r);
                    new_res->regenerate_bonds(r);
                    new_res->set_resi(
                        __original.first().first().first().size() + 1);
                    molib::Residue& add_res = chain.add(new_res);

                    // Create the relevent bond between the ligand and fragment,
                    // remove the "search atom"
                    molib::Atom& atom2 =
                        add_res.element(search_atom.atom_number());
                    molib::Atom& start2 =
                        chain.first().element(start_atom->atom_number());
                    molib::Atom& mod_atom =
                        atom2.get_bonds().front()->second_atom(atom2);
                    add_res.renumber_atoms(new_design.get_atoms().size());
                    mod_atom.connect(start2).set_bo(1);

                    new_design.set_name(
                        __original.name() + "_design_with_" + fragment.name() +
                        "_on_" + std::to_string(start_atom->atom_number()) +
                        "_to_" + std::to_string(mod_atom.atom_number()));

                    // Fix internal atom graph
                    for (size_t i = 0; i < mod_atom.size(); ++i) {
                        if (mod_atom[i].atom_number() == atom2.atom_number()) {
                            mod_atom.erase(i);
                            break;
                        }
                    }

                    mod_atom.erase_bond(atom2);
                    add_res.erase(search_atom.atom_number());
                    molib::Atom::Vec test = new_design.get_atoms();
                    molib::Atom::Set set_of_new_atoms(test.begin(), test.end());
                    if (!__existing.is_seed_unique(set_of_new_atoms)) {
                        log_note << new_design.name()
                                 << " already has been designed!" << endl;
                        __designs.erase(__designs.size() - 1);
                        continue;
                    }

                    __existing.get_seed_id(set_of_new_atoms);

                    log_note << "Created: " << new_design.name() << endl;
                    log_note << "Lipsinki is "
                             << std::get<0>(original_lipinski) +
                                    std::get<0>(frag_lipinski)
                             << " "
                             << std::get<1>(original_lipinski) +
                                    std::get<1>(frag_lipinski)
                             << " "
                             << std::get<2>(original_lipinski) +
                                    std::get<2>(frag_lipinski)
                             << endl;
                }
            }
        }
    }
}

molib::Molecules Design::functionalize_hydrogens_with_single_atoms(
    const molib::Molecule& original, const std::string& atom_type) {
    molib::Molecules designs;

    for (const auto start_atom : original.get_atoms()) {
        // Make sure the ligand atom has an open valency
        if (start_atom->get_num_hydrogens() == 0 ||
            start_atom->element() == molib::Element::H) {
            continue;
        }

        size_t count = 0;
        for (const auto bond : start_atom->get_bonds()) {
            molib::Atom& other_orig = bond->second_atom(*start_atom);

            if (other_orig.element() != molib::Element::H) continue;

            // Copy the new molecule into the returnable object
            designs.add(new molib::Molecule(original));
            designs.last().set_name(original.name() + "_added_" + atom_type +
                                    "_on_" +
                                    std::to_string(start_atom->atom_number()) +
                                    "_count_" + std::to_string(++count));

            molib::Atom& mod_atom =
                designs.last().last().last().last().last().atom(
                    other_orig.atom_number());

            mod_atom.set_element(atom_type);
            mod_atom.set_atom_name(atom_type + "_d");
            if (atom_type == "C" || atom_type == "O" || atom_type == "N" ||
                atom_type == "S") {
                mod_atom.set_idatm_type(atom_type + "3");
            } else {
                mod_atom.set_idatm_type(atom_type + "3");
            }

            geometry::Coordinate crd = (mod_atom.crd() - start_atom->crd());
            crd.normalize();
            mod_atom.set_crd(crd * 1.54 + start_atom->crd());

            // Don't generate a bunch of junk!
            if (start_atom->get_num_hydrogens() == 3) {
                break;
            }
        }
    }

    return designs;
}

void Design::functionalize_extremes_with_single_atoms(
    const std::vector<std::string>& idatms) {
    if (idatms.empty())
        throw Error("No atom types given for addition to the molecule");

    for (auto& atom_type : idatms) {
        for (auto& start_atom : __original.get_atoms()) {
            // Make sure the ligand atom has an open valency
            if (start_atom->element() == molib::Element::H ||
                start_atom->get_bonds().size() -
                        start_atom->get_num_hydrogens() !=
                    1) {
                continue;
            }

            // Copy the new molecule into the returnable object
            __designs.add(new molib::Molecule(__original));
            __designs.last().set_name(
                __original.name() + "_changed_" +
                std::to_string(start_atom->atom_number()) + "_" + atom_type);

            // Remove all hydrogens from the original ligand (they are not
            // needed anymore)
            molib::BondOrder::compute_bond_order(__designs.last().get_atoms());
            __designs.last().first().first().first().first().erase_hydrogen();

            molib::Residue& mod_residue =
                __designs.last().first().first().first().first();

            molib::Atom& atom_to_change =
                mod_residue.element(start_atom->atom_number());
            atom_to_change.set_atom_name(atom_type);
            atom_to_change.set_element(molib::Element(atom_type));
            atom_type == "C" || atom_type == "O" || atom_type == "N" ||
                    atom_type == "S"
                ? atom_to_change.set_idatm_type(atom_type + "3")
                : atom_to_change.set_idatm_type(atom_type);
            atom_to_change.set_gaff_type("???");
        }
    }
}

void design::Design::change_original_name(const std::string& name) {
    __original.set_name(name);
}
}
}
