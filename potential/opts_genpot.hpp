#ifndef OPTS_GENPOT_H
#define OPTS_GENPOT_H
#include <tclap/CmdLine.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "helper/error.hpp"
#include <thread>
using namespace std;

class CmdLnOpts {
	
	string __ref_state;

	string __rad_or_raw;
	int __dist_cutoff;
	string __distributions_file;
	double __step_non_bond;
	double __scale_non_bond;
	string __potential_file;

	bool __quiet;
	string __program_name;
	string __version;
public:
	CmdLnOpts() : __quiet(false), __version("1.00") {}
	void init(int argc, char* argv[]) {
		__program_name = argv[0];
		try {			
			TCLAP::CmdLine cmd("Command description message", ' ', __version);
			TCLAP::SwitchArg quietSwitch("q","quiet","Quiet mode (default is verbose)", cmd, false);

			vector<string> allowedRef{"mean","cumulative"};
			TCLAP::ValuesConstraint<string> allowedValsRef( allowedRef );
			TCLAP::ValueArg<string> ref_stateArg("","ref", "Normalization method for the reference state ('mean' (default) is averaged over all atom type pairs, whereas 'cumulative' is a summation for atom type pairs)",false,"mean",&allowedValsRef, cmd);

			vector<string> allowedRad{"radial","normalized_frequency"};
			TCLAP::ValuesConstraint<string> allowedValsRad( allowedRad );
			TCLAP::ValueArg<string> rad_or_rawArg("","func","Function for calculating scores 'radial' (default) or \
				'normalized_frequency'",false,"radial",&allowedValsRad, cmd);
			vector<int> allowedDistCutoff{4,5,6,7,8,9,10,11,12,13,14,15};
			TCLAP::ValuesConstraint<int> allowedValsDistCutoff( allowedDistCutoff );
			
			TCLAP::ValueArg<int> dist_cutoffArg("","cutoff","Cutoff length (default is 15)",false,15,&allowedValsDistCutoff, cmd);

			TCLAP::ValueArg<string> distributions_fileArg("","dist","Select one of the interatomic distance distribution \
				file(s) provided with this script",false,"data/csd_complete_distance_distributions.txt","string", cmd);
			TCLAP::ValueArg<double> step_non_bondArg("","step","Step for spline generation of non-bonded knowledge-based \
				potential [0.0-1.0] (default is 0.01)",false,0.01,"double", cmd);
			TCLAP::ValueArg<double> scale_non_bondArg("","scale","Scale non-bonded forces and energy for knowledge-based \
				potential [0.0-1.0] (default is 1.0)",false,1.0,"double", cmd);

			TCLAP::ValueArg<string> potential_fileArg("","potential_file","Output file for potentials and derivatives",false,"potentials.txt","string", cmd);

			cmd.parse( argc, argv );
			__quiet = quietSwitch.getValue();

			__ref_state = ref_stateArg.getValue();
			__rad_or_raw = rad_or_rawArg.getValue();
			__dist_cutoff = dist_cutoffArg.getValue();
			__distributions_file = distributions_fileArg.getValue();
			__step_non_bond = step_non_bondArg.getValue();
			__scale_non_bond = scale_non_bondArg.getValue();
			__potential_file = potential_fileArg.getValue();

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
	bool quiet() const { return __quiet; }

	string ref_state() const { return __ref_state; }
	string rad_or_raw() const { return __rad_or_raw; }
	int dist_cutoff() const { return __dist_cutoff; }
	string distributions_file() const { return __distributions_file; }
	double step_non_bond() const { return __step_non_bond; }
	double scale_non_bond() const { return __scale_non_bond; }
	string potential_file() const { return __potential_file; }

	friend ostream& operator<< (ostream& stream, const CmdLnOpts &cmdl) {
		unsigned int n = thread::hardware_concurrency();
		stream << endl << "Detected support for " << n << " concurrent threads." << endl;
		return stream;
	}

};
extern CmdLnOpts cmdl;
#endif
