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
	string __query_file;
	string __infile;
	string __outfile;
	string __query_chain_id;
	//~ string __qpdb_file;
	//~ string __qcid;
	string __pdb_dirname;
	string __lig_code;
	string __bsite_file;
	//~ string __nrpdb_file;
	//~ string __bio;
	//~ string __models;
	//~ bool __hydrogens;
	//~ bool __neighb;
	//~ bool __cluster;
	//~ bool __noalch;
	//~ bool __ralch;
	bool __quiet;
	//~ bool __bsite;
	//~ bool __geo;
	//~ bool __rnolig;
	//~ bool __nrpdb;
	string __program_name;
	string __version;
public:
	CmdLnOpts() : __quiet(false), __version("3.00") {}
	void init(int argc, char* argv[]) {
		__program_name = argv[0];
		try {			
			TCLAP::CmdLine cmd("Command description message", ' ', __version);
			TCLAP::SwitchArg quietSwitch("q","quiet","Quiet mode (default is verbose)", cmd, false);
			//~ TCLAP::SwitchArg bsiteSwitch("","bsite","Output binding site residues too", cmd, false);
			TCLAP::ValueArg<string> bsite_fileArg("","bsitefile","Output binding site residues to this file",true,"","string", cmd);
			//~ TCLAP::ValueArg<string> nrpdb_fileArg("","nrpdbfile","Read binding site residues from this file",false,"bsite.pdb","string", cmd);
			TCLAP::ValueArg<string> outfileArg("","outfile","Output ligand to this file",true,"ligand.pdb","string", cmd);
			TCLAP::ValueArg<string> lig_codeArg("","ligcode","Ligand codes to extract", true, "", "string", cmd);
			TCLAP::ValueArg<string> pdb_dirnameArg("","pdbdir","Directory containing PDB files (default is .)",false,".","string", cmd);
			TCLAP::ValueArg<string> infileArg("","infile","JSON or NOSQL filename with ProBiS alignments",true,"","string", cmd);
			TCLAP::ValueArg<string> query_fileArg("","query","Query PDB filename",true,"","string", cmd);
			TCLAP::ValueArg<string> query_chain_idArg("","qcid","Query chain id",true,"","string", cmd);
			//~ Parse the argv array.
			cmd.parse(argc, argv);
			__outfile = outfileArg.getValue();
			__infile = infileArg.getValue();
			__query_file = query_fileArg.getValue();
			__query_chain_id = query_chain_idArg.getValue();
			__bsite_file = bsite_fileArg.getValue();
			//~ __nrpdb_file = nrpdb_fileArg.getValue();
			__lig_code = lig_codeArg.getValue();
			__pdb_dirname = pdb_dirnameArg.getValue();
			//~ __bsite = bsiteSwitch.getValue();
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
	string outfile() const {return __outfile; }
	string infile() const {return __infile; }
	string query_file() const {return __query_file; }
	string query_chain_id() const {return __query_chain_id; }
	string bsite_file() const {return __bsite_file; }
	//~ string nrpdb_file() const {return __nrpdb_file; }
	string lig_code() const { return __lig_code; }
	string pdb_dirname() const { return __pdb_dirname; }
	//~ bool bsite() const { return __bsite; }
	bool quiet() const { return __quiet; }
	string version() const { return __version; }
	string program_name() const { return __program_name; }
};
extern CmdLnOpts cmdl;
#endif
