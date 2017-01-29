#include <iostream>
#include "program/cmdlnopts.hpp"

////////////////// BINDING SITE DETECTION USING PROBIS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Program::CmdLnOpts *cmdlopts = new Program::CmdLnOpts;

                cmdlopts->init(argc, argv);
                cout << *cmdlopts << endl;
                
        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
