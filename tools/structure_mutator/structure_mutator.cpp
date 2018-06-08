#include "candock/design/design.hpp"
#include <boost/program_options.hpp>
#include "candock/fileout/fileout.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"
#include "candock/parser/fileparser.hpp"
#include "candock/fragmenter/unique.hpp"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                std::vector<std::string> inputs;
                std::vector<std::string> atom_types;

                boost::program_options::variables_map vm;
                po::options_description cmdln_options;

                po::options_description generic;
                generic.add_options()
                ("help,h", "Show this help")
                ("input,f", po::value<std::vector<std::string>>(&inputs), "Input file(s)")
                ("a,add_single_atoms",    po::value< std::vector<std::string>>(&atom_types),
                 "Change hydrogens to given atoms. Multiple atoms can be given.")
                ("concatinate,c","Concatinate multiple binned files together");

                cmdln_options.add(generic);

                po::positional_options_description p;
                p.add("input", -1);

                po::store (po::command_line_parser (argc, argv).options(cmdln_options).positional(p).run(), vm);
                po::notify (vm);

                if (vm.count ("help") || (vm.count("extract_only") && vm.count("bin_only"))) {
                        std::cout << cmdln_options << std::endl;
                        return 0;
                }

                for (auto input_file : inputs) {
                        Parser::FileParser input(
                            input_file,
                            Parser::pdb_read_options::all_models | Parser::pdb_read_options::hydrogens);

                        Molib::Molecules mols;

                        while (input.parse_molecule(mols)) {

                                mols.compute_idatm_type();
                                mols.compute_bond_order();

                                for ( const auto& mol : mols ) {
                                        Molib::Molecules mods =  design::Design::functionalize_hydrogens_with_single_atoms(mol,atom_types[0]);

                                        for (const auto& mol2 : mods)
                                                Fileout::print_mol2(std::cout, mol2);
                                }

                                mols.clear();
                        }
                }
                

        } catch (po::error& e) {
                std::cerr << "error: " << e.what() << std::endl;
                std::cerr << "Please see the help (-h) for more information" << std::endl;
        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}
