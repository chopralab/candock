#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"
#include "modeler/systemtopology.hpp"

////////////////// LINKING OF FRAGMENTS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Program::CmdLnOpts *cmdlopts = new Program::CmdLnOpts;
                cmdlopts->init(argc, argv, Program::CmdLnOpts::STARTING |
                                           Program::CmdLnOpts::FORCE_FIELD |
                                           Program::CmdLnOpts::SCORING |
                                           Program::CmdLnOpts::LINKING );

                Benchmark main_timer;
                main_timer.display_time("started");
                cout << *cmdlopts << endl;

                help::Options::set_options(cmdlopts);

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                Program::Target targets (cmdl.get_string_option("receptor"));
                targets.find_centroids();
                targets.dock_fragments(ligand_fragmenter);
                OMMIface::SystemTopology::loadPlugins();
                targets.link_fragments();

                main_timer.display_time("Finished");

        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
