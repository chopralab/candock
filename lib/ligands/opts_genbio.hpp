#ifndef OPTS_GENBIO_H
#define OPTS_GENBIO_H
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
	string __biofile;
	string __mols_name;
	string __qpdb_file;
	string __qcid;
	string __pdb_dirname;
	string __bio;
	string __models;
	bool __hydrogens;
	bool __neighb;
	bool __noalch;
	bool __ralch;
	bool __quiet;
	bool __rnolig;
	bool __rasym;
	string __program_name;
	string __version;
public:
	CmdLnOpts() : __quiet(false), __version("3.00") {}
	void init(int argc, char* argv[]) {
		__program_name = argv[0];
		try {			
			TCLAP::CmdLine cmd("Command description message", ' ', __version);
			TCLAP::SwitchArg quietSwitch("q","quiet","Quiet mode (default is verbose)", cmd, false);
			TCLAP::ValueArg<string> biofileArg("","biofile","Output biological units to this pdb file",false,"","string", cmd);
			TCLAP::ValueArg<string> mols_nameArg("","mols_name","Stem name of representative",false,"","string", cmd);
			vector<string> allowedModels{"all","first"};
			TCLAP::ValuesConstraint<string> allowedValsModels( allowedModels );
			TCLAP::ValueArg<string> modelsArg("","models","Read models (default is all)", false, "all", &allowedValsModels, cmd);
			vector<string> allowedBio{"all","first", "none"};
			TCLAP::ValuesConstraint<string> allowedValsBio( allowedBio );
			TCLAP::ValueArg<string> bioArg("","bio","Generate biological assemblies (default is none)", false, "none", &allowedValsBio, cmd);
			TCLAP::ValueArg<string> qcidArg("","qcid","Query chain id(s) (default is A)", false, "A", "string", cmd);
			TCLAP::ValueArg<string> qpdb_fileArg("","qpdb","Query pdb file", false, "", "string", cmd);
			TCLAP::SwitchArg hydrogensSwitch("","hydro","Read hydrogens", cmd, false);
			TCLAP::SwitchArg neighbSwitch("","neighb","Find and output (in REMARK) neighbor residues of ligands", cmd, false);
			TCLAP::SwitchArg rnoligSwitch("","rnolig","Remove all aligned chains that are not ligands (within 4.0 A) of the query protein", cmd, false);
			TCLAP::SwitchArg rasymSwitch("","rasym","Remove asymmetric unit when biological assembly exists", cmd, false);
			TCLAP::SwitchArg ralchSwitch("","ralch","Remove assemblies that don't have aligned chains", cmd, false);
			TCLAP::SwitchArg noalchSwitch("","noalch","Remove aligned protein chains (query chain and aligned chains) from the first model of each assembly", cmd, false);
			TCLAP::ValueArg<string> pdb_dirnameArg("","pdbdir","Directory containing PDB files (default is .)",false,".","string", cmd);
			TCLAP::ValueArg<string> infileArg("","infile","JSON or NOSQL filename with ProBiS alignments",false,"","string", cmd);
			//~ Parse the argv array.
			cmd.parse(argc, argv);
			__biofile = biofileArg.getValue();
			__mols_name = mols_nameArg.getValue();
			__infile = infileArg.getValue();
			__bio = bioArg.getValue();
			__models = modelsArg.getValue();
			__qcid = qcidArg.getValue();
			__qpdb_file = qpdb_fileArg.getValue();
			__pdb_dirname = pdb_dirnameArg.getValue();
			__hydrogens = hydrogensSwitch.getValue();
			__neighb = neighbSwitch.getValue();
			__noalch = noalchSwitch.getValue();
			__ralch = ralchSwitch.getValue();
			__quiet = quietSwitch.getValue();
			__rnolig = rnoligSwitch.getValue();
			__rasym = rasymSwitch.getValue();
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
	string biofile() const {return __biofile; }
	string mols_name() const {return __mols_name; }
	string infile() const {return __infile; }
	string qpdb_file() const { return __qpdb_file; }
	string qcid() const { return __qcid; }
	string bio() const { return __bio; }
	string models() const { return __models; }
	string pdb_dirname() const { return __pdb_dirname; }
	bool hydrogens() const { return __hydrogens; }
	bool neighb() const { return __neighb; }
	bool noalch() const { return __noalch; }
	bool ralch() const { return __ralch; }
	bool quiet() const { return __quiet; }
	bool rnolig() const { return __rnolig; }
	bool rasym() const { return __rasym; }
	string version() const { return __version; }
	string program_name() const { return __program_name; }
};
extern CmdLnOpts cmdl;
#endif
