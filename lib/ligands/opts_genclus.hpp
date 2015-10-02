#ifndef OPTS_GENCLUS_H
#define OPTS_GENCLUS_H
#include <tclap/CmdLine.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "helper/error.hpp"
using namespace std;

class CmdLnOpts {
	string __json_file;
	string __json_with_ligs_file;
	string __bio_dir;
	string __names_dir;
	bool __neighb;
	double __probis_clus_rad;
	int __probis_min_pts;
	double __probis_min_z_score;

	bool __quiet;
	string __program_name;
	string __version;
public:
	CmdLnOpts() : __quiet(false), __version("3.00") {}
	void init(int argc, char* argv[]) {
		__program_name = argv[0];
		try {			
			TCLAP::CmdLine cmd("Command description message", ' ', __version);
			TCLAP::SwitchArg quietSwitch("q","quiet","Quiet mode (default is verbose)", cmd, false);
			TCLAP::ValueArg<string> json_fileArg("","json","Json-formatted ProBiS alignments input file (default is probis.json)",true,"probis.json","string", cmd);
			TCLAP::ValueArg<string> json_with_ligs_fileArg("","jsonwl","Json-formatted ProBiS alignments with transposed ligands output file (default is probis_with_ligands.json",false,"probis_with_ligands.json","string", cmd);
			TCLAP::ValueArg<string> bio_dirArg("","bio","Directory with ProBiS-ligands bio database (default is data/probis_ligands/bio)",false,"data/probis_ligands/bio","string", cmd);
			TCLAP::ValueArg<string> names_dirArg("","names","Directory with ligand names (default is data/probis_ligands/names)",false,"data/probis_ligands/names","string", cmd);
			TCLAP::SwitchArg neighbSwitch("","neighb","Allow only ligands that are in the similar regions according to REMARKs (not enabled by default)", cmd, false);
			TCLAP::ValueArg<double> probis_clus_radArg("","probis_clus_rad","Cluster radius for predicted ligands by probis (default is 2.0)",false,2.0,"double", cmd);
			TCLAP::ValueArg<int> probis_min_ptsArg("","probis_min_pts","The minimum number of points (for predicted ligands) required to form a cluster (default is 10)",false,10,"int", cmd);
			TCLAP::ValueArg<double> probis_min_z_scoreArg("","probis_min_z_score","Minimium z-score of ligands to be considered in clustering (default is 2.0)",false,2.0,"double", cmd);
			//~ Parse the argv array.
			cmd.parse(argc, argv);
			__json_file = json_fileArg.getValue();
			__json_with_ligs_file = json_with_ligs_fileArg.getValue();
			__bio_dir = bio_dirArg.getValue();
			__names_dir = names_dirArg.getValue();
			__neighb = neighbSwitch.getValue();
			__probis_clus_rad = probis_clus_radArg.getValue();
			__probis_min_pts = probis_min_ptsArg.getValue();
			__probis_min_z_score = probis_min_z_scoreArg.getValue();
		} 
		catch (TCLAP::ArgException &e) { 
			cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
			throw Error("die: arguments error\n");
		}		
	}
	void display_time(string what) {
		cout << what << " on " << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "\n";
		cout << "running " << __program_name << " version " << __version << " on hostname " << boost::asio::ip::host_name() << "\n";
	}
	// interface
	string json_file() const { return __json_file; }
	string json_with_ligs_file() const { return __json_with_ligs_file; }
	string bio_dir() const { return __bio_dir; }
	string names_dir() const { return __names_dir; }
	bool neighb() const { return __neighb; }
	double probis_clus_rad() const { return __probis_clus_rad; }
	int probis_min_pts() const { return __probis_min_pts; }
	double probis_min_z_score() const { return __probis_min_z_score; }

	string version() const { return __version; }
	string program_name() const { return __program_name; }
};
extern CmdLnOpts cmdl;
#endif
