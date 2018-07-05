#include "candock/probis/probis.hpp"
#include "const.h"
#include "args.h"
#include "states.h"
#include "candock/helper/error.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/benchmark.hpp"

using namespace std;

namespace probis {
	void compare_against_bslib(const string &receptor_file, const string &surface_file,
		const string &receptor_chain_id, const string &bslib_file, const int ncpu, 
		const string &nosql_file, const string &json_file) {
		try {
			Benchmark bench;
			cout << "Starting ProBiS for binding sites prediction" << endl;
		    Args args;
		    _longnames = true;
		    _local = true;
		    //~ args.fill_args(argc, argv);
			NCPU = ncpu;
			SRF_FILE = surface_file;
			PROTEIN1 = receptor_file;
			CHAIN1 = receptor_chain_id;
			//~ strcpy(SURF_FILE, bslib_file.c_str());
			SURF_FILE = bslib_file.c_str();
			NOSQL_FILE = nosql_file;
			JSON_FILE = json_file;
		    args.print_state();
		    args.print_modifiers();
		    args.print_constants();
			state3(&args);
			_srf = true;
			PROTEIN1 = SRF_FILE;
			state4(&args);
			PROTEIN1 = receptor_file;
			state7(&args);
			cout << "Binding sites prediction took " 
				<< bench.seconds_from_start() << " wallclock seconds\n";

		}  // .. konec try bloka
		catch (Err e) {
			cout << e.what() << endl;
			throw Error("die : something went wrong in probis ...");
		}
	}
	
};
