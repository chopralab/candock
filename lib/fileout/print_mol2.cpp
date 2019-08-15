/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "fileout.hpp"

namespace Fileout {
    
        const map<const string, const string> id_atom_to_sybyl {
                {"C3", "C.3"},
                {"C2", "C.2"},
                {"Cac","C.2"},
                {"Car","C.ar"},
                {"C1", "C.1"},
                {"N3+","N.3"},
                {"N3", "N.3"},
                {"N2", "N.2"},
                {"N1", "N.1"},
                {"Ntr","N.2"},
                {"O3", "O.3"},
                {"O3-","O.3"},
                {"O2", "O.2"},
                {"O2-","O.2"},
                {"S3", "S.3"},
                {"Nar","N.ar"},
                {"Pac","P.3"},
                {"Pax","P.3"},
                {"P",  "P.3"},
                {"S3", "S.3"},
                {"S2", "S.2"},
                {"Npl","N.pl3"},
                {"???","ANY"},
        };

    
        void print_mol2( std::ostream &ss, const Molib::Molecule &ligand ) {
                const Molib::Atom::Vec all_atoms = ligand.get_atoms();
                const Molib::BondSet all_bonds = Molib::get_bonds_in(all_atoms);

                ss << "@<TRIPOS>MOLECULE\n";
                ss << ligand.name() << "\n";
                ss << all_atoms.size() << " " << all_bonds.size() << " 1\n";
                ss << "SMALL\nUSER_CHARGES\n";

                std::map <int, int> atom_mapper;
                ss << "@<TRIPOS>ATOM\n";

                for (size_t i = 0; i < all_atoms.size(); ++i) {
                        atom_mapper.insert( {all_atoms[i]->atom_number(), i + 1} );

                        ss << std::setw(10) << std::left << i + 1;
                        ss << std::setw(5)  << all_atoms[i]->atom_name();

                        ss << std::setw(15) << std::setprecision(5) << all_atoms[i]->crd().x();
                        ss << std::setw(15) << std::setprecision(5) << all_atoms[i]->crd().y();
                        ss << std::setw(15) << std::setprecision(5) << all_atoms[i]->crd().z();

                        if ( id_atom_to_sybyl.count( all_atoms[i]->idatm_type_unmask() ) != 0 ) {
                                ss << std::setw(5) << id_atom_to_sybyl.at(all_atoms[i]->idatm_type_unmask());
                        } else {
                                ss << std::setw(5) << all_atoms[i]->idatm_type_unmask();
                        }

                        ss << " 1 ";
                        ss << all_atoms[i]->br().resn() + std::to_string(all_atoms[i]->br().resi());
                        ss << " 0.0000\n";
                }

                ss << "@<TRIPOS>BOND\n";
                size_t counter = 0;
                for ( const Molib::Bond* b : all_bonds) {
                        ss << setw(10) << ++counter;
                        ss << setw(10) << atom_mapper.at(b->atom1().atom_number());
                        ss << setw(10) << atom_mapper.at(b->atom2().atom_number());
                        ss << " " << (b->get_bo() == 0? 1 : b->get_bo()) << "\n";
                }
                ss << "\n";
        }
}
