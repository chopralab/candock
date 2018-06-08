#include <iostream>
#include "candock/program/cmdlnopts.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/program/target.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"

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

                Score::KBFF score(cmdl.get_string_option("ref"), "complete",
                                   cmdl.get_string_option("func"), cmdl.get_int_option("cutoff"),
                                   cmdl.get_double_option("step"));

                // Pretend that we'll use all types
                set<int> all_types;
                int sz = help::idatm_mask.size();
                for (int i = 0; i < sz; ++i) {
                        all_types.insert(i);
                }

                score.define_composition(all_types, all_types)
                     .process_distributions_file(cmdl.get_string_option("dist"));

                score.compile_objective_function(cmdl.get_double_option("scale"));

                score.output_objective_function(cmdl.get_string_option("obj_dir"));

                Inout::output_file(score, cmdl.get_string_option("potential_file"));

                main_timer.display_time("Finished");
        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}
