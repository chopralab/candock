#include <iostream>
#include "program/cmdlnopts.hpp"
#include "helper/benchmark.hpp"
#include "program/target.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

using namespace std;

////////////////// GENERATE POTENTIAL FUNCTIONS ///////////////////////////

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                help::Options::set_options( new Program::CmdLnOpts( 
                    argc, argv, Program::CmdLnOpts::SCORING));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                cout << help::Options::get_options()->configuration_file() << endl;

                Score::Score score(cmdl.get_string_option("ref"), "complete",
                                   cmdl.get_string_option("func"), cmdl.get_int_option("cutoff"),
                                   cmdl.get_double_option("step"));

                score.define_composition(set<int>(), set<int>())
                     .process_distributions_file(cmdl.get_string_option("dist"))
                     .compile_objective_function(cmdl.get_double_option("scale"));

                score.output_objective_function(cmdl.get_string_option("obj_dir"));

                Inout::output_file(score, cmdl.get_string_option("potential_file"));

                main_timer.display_time("Finished");
        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}
