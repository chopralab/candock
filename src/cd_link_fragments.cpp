#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"
#include "modeler/systemtopology.hpp"
#include "version.hpp"


////////////////// LINKING OF FRAGMENTS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                help::Options::set_options( new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::FORCE_FIELD |
                                Program::CmdLnOpts::SCORING |
                                Program::CmdLnOpts::LINKING ));

                Benchmark main_timer;
                main_timer.display_time("started");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                log_note << help::Options::get_options()->configuration_file() << endl;

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                Program::Target targets (cmdl.get_string_option("receptor"));
                targets.link_fragments(ligand_fragmenter);

                main_timer.display_time("Finished");

        } catch (exception& e) {
                log_error << e.what() << endl;
                return 1;
        }
        return 0;
}
