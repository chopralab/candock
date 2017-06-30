#include "program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {
                if ( argc <= 1 ) {
                        cerr << "You MUST supply an argument!" << endl;
                        return 1;
                }

                Inout::Logger::set_all_stderr(true);

                Parser::FileParser input(argv[1], Parser::pdb_read_options::all_models);

                Molib::Molecules input_read;

                input.parse_molecule(input_read);

                input_read.compute_idatm_type()
                .compute_hydrogen()
                .compute_bond_order()
                .compute_bond_gaff_type()
                .refine_idatm_type()
                .erase_hydrogen()  // needed because refine changes connectivities
                .compute_hydrogen()   // needed because refine changes connectivities
                .compute_ring_type()
                .compute_gaff_type()
                .compute_rotatable_bonds() // relies on hydrogens being assigned
                .erase_hydrogen();

                cout << input_read;

        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}

