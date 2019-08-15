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
                .compute_ring_type()
                .compute_gaff_type()
                .compute_rotatable_bonds(); // relies on hydrogens being assigned

                cout << input_read;

        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}

