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
                __stream.seekg(pos);
        }
        
        void Parser::set_hm(unsigned int hm) {
            __hm = hm;
        }
}
