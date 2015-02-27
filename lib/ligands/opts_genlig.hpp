#ifndef CMDLNOPTS_H
#define CMDLNOPTS_H
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

	string __receptor_file;
	string __ligand_file;
	string __receptor_chain_id;

	string __json_file;
	string __bio_dir;

	string __lig_code;
	string __bsite_file;

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

			TCLAP::ValueArg<string> receptor_fileArg("","receptor","Receptor filename",true,"default","string", cmd);
			TCLAP::ValueArg<string> ligand_fileArg("","ligand","Ligand filename",true,"default","string", cmd);
			TCLAP::ValueArg<string> receptor_chain_idArg("","receptor_chain_id","Chain id(s) of the receptor (default is A)",
				false,"A","string", cmd);

			TCLAP::ValueArg<string> json_fileArg("","json","Json-formatted ProBiS alignments output file (default is probis.json",false,"probis.json","string", cmd);
			TCLAP::ValueArg<string> bio_dirArg("","bio","Directory with ProBiS-ligands bio database (default is data/probis_ligands/bio)",false,"data/probis_ligands/bio","string", cmd);

			TCLAP::ValueArg<string> lig_codeArg("","ligcode","Ligand codes to extract", true, "", "string", cmd);
			TCLAP::ValueArg<string> bsite_fileArg("","bsitefile","Output binding site residues to this file",true,"","string", cmd);

			//~ Parse the argv array.
			cmd.parse(argc, argv);

			__receptor_file = receptor_fileArg.getValue();
			__ligand_file = ligand_fileArg.getValue();
			__receptor_chain_id = receptor_chain_idArg.getValue();

			__json_file = json_fileArg.getValue();
			__bio_dir = bio_dirArg.getValue();

			__lig_code = lig_codeArg.getValue();
			__bsite_file = bsite_fileArg.getValue();

			__quiet = quietSwitch.getValue();
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
	string receptor_file() const { return __receptor_file; }
	string ligand_file() const { return __ligand_file; }
	string receptor_chain_id() const { return __receptor_chain_id; }

	string json_file() const { return __json_file; }
	string bio_dir() const { return __bio_dir; }

	string lig_code() const { return __lig_code; }
	string bsite_file() const {return __bsite_file; }

	bool quiet() const { return __quiet; }
	string version() const { return __version; }
	string program_name() const { return __program_name; }
};
extern CmdLnOpts cmdl;
#endif
