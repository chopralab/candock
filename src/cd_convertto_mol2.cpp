#include "candock/program/cmdlnopts.hpp"
#include "candock/program/target.hpp"

#include "candock/score/score.hpp"
#include "candock/modeler/forcefield.hpp"
#include "candock/modeler/modeler.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"

#include "candock/fileout/fileout.hpp"

using namespace std;
using namespace candock;
using namespace candock::Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                for (int i = 1; i < argc; ++i) {
                        parser::FileParser input(argv[i], parser::pdb_read_options::all_models);
                        molib::Molecules mols;
                        input.parse_molecule(mols);

                        for (const auto &m : mols) {
                                fileout::print_mol2(std::cout,m);
                        }
                }


        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}


