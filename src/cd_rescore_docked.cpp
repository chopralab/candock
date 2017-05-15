#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

using namespace std;
using namespace Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
        try {

                if(!drm::check_drm()) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                help::Options::set_options(new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::LIG_FRAMGENT| 
                                Program::CmdLnOpts::SCORING));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                log_note << help::Options::get_options()->configuration_file() << endl;


                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                Program::Target targets(cmdl.get_string_option("receptor"));
                targets.rescore_docked(ligand_fragmenter);

        } catch (exception& e) {
                log_error << e.what() << endl;
        }
        return 0;
}

