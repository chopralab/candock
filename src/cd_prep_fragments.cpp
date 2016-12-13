#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
	try {

		Program::CmdLnOpts cmdl;

		cmdl.init(argc, argv, Program::CmdLnOpts::STARTING |
		                      Program::CmdLnOpts::LIG_FRAMGENT );
		cout << cmdl << endl;

		cmdl.display_time("Starting");

		Program::FragmentLigands ligand_fragmenter;
		ligand_fragmenter.run_step(cmdl);

		cmdl.display_time("Finished");

	} catch (exception& e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}
