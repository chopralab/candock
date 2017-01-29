#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/target.hpp"

////////////////// BINDING SITE DETECTION USING PROBIS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Program::CmdLnOpts *cmdlopts = new Program::CmdLnOpts;
                cmdlopts->init(argc, argv, Program::CmdLnOpts::STARTING |
                                           Program::CmdLnOpts::PROBIS);

                Benchmark main_timer;
                main_timer.display_time("Started");
                cout << *cmdlopts << endl;

                help::Options::set_options(cmdlopts);

                Program::Target targets (cmdl.get_string_option("receptor"));
                targets.find_centroids();

                main_timer.display_time("Finished");


        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
