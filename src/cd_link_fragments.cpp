#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"

////////////////// LINKING OF FRAGMENTS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
	try {

		Program::CmdLnOpts cmdl;

		cmdl.init(argc, argv, Program::CmdLnOpts::STARTING |
		                      Program::CmdLnOpts::FORCE_FIELD |
		                      Program::CmdLnOpts::SCORING |
		                      Program::CmdLnOpts::LINKING );
		cout << cmdl << endl;

		cmdl.display_time("Starting");

		Program::FragmentLigands ligand_fragmenter;
		ligand_fragmenter.run_step(cmdl);

		Program::Target targets (cmdl.receptor_file());
		targets.find_centroids(cmdl);
		targets.dock_fragments(ligand_fragmenter, cmdl);
		targets.link_fragments(cmdl);

		cmdl.display_time("Finished");

	} catch (exception& e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}
