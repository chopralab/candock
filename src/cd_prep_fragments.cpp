#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Program::CmdLnOpts *cmdlopts = new Program::CmdLnOpts;
                cmdlopts->init(argc, argv, Program::CmdLnOpts::STARTING |
                                           Program::CmdLnOpts::LIG_FRAMGENT );

                cout << *cmdlopts << endl;

                Benchmark main_timer;
                main_timer.display_time("started");
                cout << *cmdlopts << endl;

                help::Options::set_options(cmdlopts);

                main_timer.display_time("Starting");

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                main_timer.display_time("Finished");

        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
