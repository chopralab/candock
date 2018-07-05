#include "candock/fileout/fileout.hpp"

namespace candock{
namespace fileout {
    
        const map<const string, const string> id_atom_to_sybyl {
                {"C3", "C.3"},
                {"C2", "C.2"},
                {"Cac","C.2"},
                {"Car","C.ar"},
                {"C1", "C.1"},
                {"HC", "H"},
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

        void print_mol2( std::ostream &ss, const molib::Molecule &ligand ) {
                const molib::Atom::Vec all_atoms = ligand.get_atoms();
                const molib::BondSet all_bonds = molib::get_bonds_in(all_atoms);

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
                for ( const molib::Bond* b : all_bonds) {
                        ss << setw(10) << ++counter;
                        ss << setw(10) << atom_mapper.at(b->atom1().atom_number());
                        ss << setw(10) << atom_mapper.at(b->atom2().atom_number());
                        ss << " " << (b->get_bo() == 0? 1 : b->get_bo()) << "\n";
                }
                ss << "\n";
        }
}
}
