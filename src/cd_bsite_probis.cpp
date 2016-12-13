#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/target.hpp"

////////////////// BINDING SITE DETECTION USING PROBIS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
	try {

		Program::CmdLnOpts cmdl;

		cmdl.init(argc, argv, Program::CmdLnOpts::STARTING |
		                      Program::CmdLnOpts::PROBIS);
		cout << cmdl << endl;

		cmdl.display_time("Starting");

		Program::Target targets (cmdl.receptor_file());
		targets.find_centroids(cmdl);

		cmdl.display_time("Finished");

	} catch (exception& e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}
