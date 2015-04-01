#include <iostream>
#include <exception>
#include <typeinfo>
#include "jsonreader.hpp"
#include "nosqlreader.hpp"
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"
#include "opts_genclus.hpp"
#include "geom3d/matrix.hpp"
#include "helper/error.hpp"
#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "pdbreader/grid.hpp"
#include "common.hpp"
#include "genclus.hpp"
using namespace std;

CmdLnOpts cmdl;

// g++ -static -std=c++0x gclus.cpp -o gclus -L. -lligands -L../pdbreader -lpdb -L../geom3d -lgeom3d -L../helper -lhelper -L../fragmenter -lfragmenter -L../pdbreader -lpdb -I.. -I../tclap-1.2.1/include -I../jsoncpp-src-0.6.0-rc2/include -I../jsoncpp-src-0.6.0-rc2/src -lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time -lgsl -lgslcblas -lm
int main(int argc, char *argv[]) {
	try {
		cmdl.init(argc, argv);
		genclus::generate_clusters_of_ligands(cmdl.json_file(), cmdl.json_with_ligs_file(),
			cmdl.geo_dir(), cmdl.names_dir(), cmdl.neighb(), cmdl.probis_clus_rad(),
			cmdl.probis_min_pts(), cmdl.probis_min_z_score(), true);
		cout << "end" << endl;
	}
	catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
