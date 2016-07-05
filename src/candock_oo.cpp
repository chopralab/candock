#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "opts_candock.hpp"
#include "findcentroidsstep.hpp"
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecules.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		/* Create empty output files
		 * 
		 */
		inout::output_file("", cmdl.gridpdb_hcp_file()); // gridpoints for all binding sites
		inout::output_file("", cmdl.prep_file()); // output prepared ligands
		inout::output_file("", cmdl.docked_file()); // output docked molecule conformations
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments
		
		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

		Program::FindCentroidsStep find_centroids(receptors[0]);
		
	} catch ( exception& e) {
		cerr << e.what() << endl;
	}
}
