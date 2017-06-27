#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"
#include "version.hpp"
#include "drm/drm.hpp"


////////////////// LINKING OF FRAGMENTS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                if(Version::drm_active() && !drm::check_drm()) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                help::Options::set_options( new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::FORCE_FIELD |
                                Program::CmdLnOpts::SCORING |
                                Program::CmdLnOpts::LINKING ));

                Benchmark main_timer;
                main_timer.display_time("started");

                std::cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                cout << help::Options::get_options()->configuration_file() << endl;

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                Program::Target targets (cmdl.get_string_option("receptor"));
                targets.link_fragments(ligand_fragmenter);

                main_timer.display_time("Finished");

        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
