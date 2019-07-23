/* This is structure_mutator.cpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

#include <iostream>
#include <boost/program_options.hpp>

#include "candock/design/design.hpp"

#include "statchem/fileio/fileout.hpp"
#include "statchem/fragmenter/unique.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/helper/logger.hpp"
#include "version.hpp"

using namespace std;
using namespace candock;
using namespace statchem;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    try {
        Logger::set_all_stderr(true);

        std::vector<std::string> inputs;
        std::vector<std::string> atom_types;

        boost::program_options::variables_map vm;
        po::options_description cmdln_options;

        po::options_description generic;
        generic.add_options()("help,h", "Show this help")(
            "input,f", po::value<std::vector<std::string>>(&inputs),
            "Input file(s)")(
            "a,add_single_atoms",
            po::value<std::vector<std::string>>(&atom_types),
            "Change hydrogens to given atoms. Multiple atoms can be given.")(
            "concatinate,c", "Concatinate multiple binned files together");

        cmdln_options.add(generic);

        po::positional_options_description p;
        p.add("input", -1);

        po::store(po::command_line_parser(argc, argv)
                      .options(cmdln_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);

        if (vm.count("help") ||
            (vm.count("extract_only") && vm.count("bin_only"))) {
            std::cout << cmdln_options << std::endl;
            return 0;
        }

        for (auto input_file : inputs) {
            parser::FileParser input(input_file,
                                     parser::pdb_read_options::all_models |
                                         parser::pdb_read_options::hydrogens);

            molib::Molecules mols;

            while (input.parse_molecule(mols)) {
                mols.compute_idatm_type();
                mols.compute_bond_order();

                for (const auto& mol : mols) {
                    molib::Molecules mods = design::Design::
                        functionalize_hydrogens_with_single_atoms(
                            mol, atom_types[0]);

                    for (const auto& mol2 : mods)
                        fileio::print_mol2(std::cout, mol2);
                }

                mols.clear();
            }
        }

    } catch (po::error& e) {
        std::cerr << "error: " << e.what() << std::endl;
        std::cerr << "Please see the help (-h) for more information"
                  << std::endl;
    } catch (exception& e) {
        cerr << e.what() << endl;
    }
    return 0;
}
