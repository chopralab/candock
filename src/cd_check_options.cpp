#include <iostream>
#include "program/cmdlnopts.hpp"

////////////////// BINDING SITE DETECTION USING PROBIS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
	try {

		Program::CmdLnOpts cmdl;

		cmdl.init(argc, argv);
		cout << cmdl << endl;
                
	} catch (exception& e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}
