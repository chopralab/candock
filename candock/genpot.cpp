#include <set>
#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "opts_genpot.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "helper/error.hpp"
#include "score/score.hpp"
using namespace std;

CmdLnOpts cmdl;

////////////////// GENERATE POTENTIAL FUNCTIONS ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		Molib::Score score(cmdl.ref_state(), "complete", 
			cmdl.rad_or_raw(), cmdl.dist_cutoff(), cmdl.step_non_bond());

		score.define_composition(set<int>(), set<int>())
			.process_distributions_file(cmdl.distributions_file())
			.compile_objective_function()
			.output_objective_function(cmdl.obj_dir());
		
		inout::output_file(score, cmdl.potential_file());
		
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
