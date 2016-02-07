#include "opts_genbio.hpp"
#include "genbio.hpp"
using namespace std;

// g++ -static -std=c++0x gbio.cpp -o gbio -L. -lligands -L../pdbreader -lpdb -L../geom3d -lgeom3d -L../helper -lhelper -L../fragmenter -lfragmenter -L../pdbreader -lpdb -I.. -I../tclap-1.2.1/include -I../jsoncpp-src-0.6.0-rc2/include -I../jsoncpp-src-0.6.0-rc2/src -lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time -lgsl -lgslcblas -lm
CmdLnOpts cmdl;

int main(int argc, char *argv[]) {
	try {
		cmdl.init(argc, argv);
		genbio::generate_biological_assemblies(cmdl.models(), cmdl.hydrogens(),
			cmdl.pdb_dirname(), cmdl.qpdb_file(), cmdl.qcid(), cmdl.neighb(),
			cmdl.rnolig(), cmdl.bio(), cmdl.ralch(), cmdl.infile(), cmdl.biofile(),
			cmdl.noalch(), cmdl.mols_name(), cmdl.rasym());
	}
	catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
