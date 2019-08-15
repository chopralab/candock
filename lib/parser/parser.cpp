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

#include "parser.hpp"

using namespace Molib;

namespace Parser {

        void Parser::__generate_molecule(Molecules &mols, bool &found_molecule, const std::string &name) {
                // if there were no REMARK or BIOMOLECULE...
                if(!found_molecule) {

                        Molecule* mol = nullptr;

                        if (name.empty()) {
                                log_warning << "Warning: unlabled molecule. Naming: unlabeled_" << mols.size() << endl;
                                mol = new Molecule(std::string("unlabled_") + std::to_string(mols.size()) );
                        } else {
                                mol = new Molecule(name);
                        }

                        mols.add(mol);
                        found_molecule = true;
                }
        }

        void Parser::__generate_assembly(Molecules &mols, bool &found_assembly, int assembly_number, const std::string &name) {
                // if there were no REMARK or BIOMOLECULE...
                if(!found_assembly) {
                        mols.last().add(new Assembly(assembly_number, name));
                        found_assembly = true;
                }
        }

        void Parser::__generate_model(Molecules &mols, bool &found_model, int model_number) {
                // if there were no REMARK or BIOMOLECULE...
                if(!found_model) {
                        mols.last().last().add(new Model(model_number));
                        found_model = true;
                }
        }

        void Parser::set_pos(std::streampos pos) {
                __pos = pos;
        }
        
        void Parser::set_hm(unsigned int hm) {
            __hm = hm;
        }
}
