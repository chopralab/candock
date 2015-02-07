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
	string __infile;
	string __pdb_dirname;
	string __names_dir;
	bool __quiet;
	bool __neighb;
	string __program_name;
	string __version;
public:
	CmdLnOpts() : __quiet(false), __version("3.00") {}
	void init(int argc, char* argv[]) {
		__program_name = argv[0];
		try {			
			TCLAP::CmdLine cmd("Command description message", ' ', __version);
			TCLAP::SwitchArg quietSwitch("q","quiet","Quiet mode (default is verbose)", cmd, false);
			TCLAP::SwitchArg neighbSwitch("","neighb","Allow only ligands that are in the similar regions according to REMARKs", cmd, false);
			TCLAP::ValueArg<string> pdb_dirnameArg("","pdbdir","Directory containing PDB files (default is .)",false,".","string", cmd);
			TCLAP::ValueArg<string> names_dirArg("","lndir","Directory with names of ligands (default is .)",false,".","string", cmd);
			TCLAP::ValueArg<string> infileArg("","infile","JSON or NOSQL filename with ProBiS alignments",true,"","string", cmd);
			//~ Parse the argv array.
			cmd.parse(argc, argv);
			__infile = infileArg.getValue();
			__pdb_dirname = pdb_dirnameArg.getValue();
			__names_dir = names_dirArg.getValue();
			__quiet = quietSwitch.getValue();
			__neighb = neighbSwitch.getValue();
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
	string infile() const {return __infile; }
	string pdb_dirname() const { return __pdb_dirname; }
	string names_dir() const { return __names_dir; }
	bool quiet() const { return __quiet; }
	bool neighb() const { return __neighb; }
	string version() const { return __version; }
	string program_name() const { return __program_name; }
};
extern CmdLnOpts cmdl;
#endif
