#include "opts_genlig.hpp"
#include "genlig.hpp"
using namespace std;

// g++ -static -std=c++0x glig.cpp -o glig -L. -lligands -L../pdbreader -lpdb -L../geom3d -lgeom3d -L../helper -lhelper -L../fragmenter -lfragmenter -L../pdbreader -lpdb -I.. -I../tclap-1.2.1/include -I../jsoncpp-src-0.6.0-rc2/include -I../jsoncpp-src-0.6.0-rc2/src -lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time -lgsl -lgslcblas -lm
CmdLnOpts cmdl;

int main(int argc, char *argv[]) {
	try {
		cmdl.init(argc, argv);
		genlig::generate_ligands(cmdl.receptor_file(), cmdl.receptor_chain_id(),
		cmdl.json_file(), cmdl.bio_dir(), cmdl.lig_code(),
		cmdl.ligand_file(), cmdl.bsite_file());
	}
	catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
