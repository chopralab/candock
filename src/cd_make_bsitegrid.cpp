#include "candock/program/fragmentligands.hpp"
#include "candock/program/findcentroids.hpp"
#include "candock/program/target.hpp"
#include "candock/drm/drm.hpp"
#include "version.hpp"

////////////////// PRODUCE A GRIDHCP FILE IN PDB FORMAT ///////////////////////////

using namespace std;
using namespace candock;

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                help::Options::set_options( new Program::CmdLnOpts (
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::SCORING |
                                Program::CmdLnOpts::FRAG_DOCKING ));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                cout << help::Options::get_options()->configuration_file() << endl;

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                Program::Target targets (cmdl.get_string_option("receptor"));
                targets.make_gridhcp(ligand_fragmenter);

        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
