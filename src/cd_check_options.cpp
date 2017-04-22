#include <iostream>
#include "program/cmdlnopts.hpp"
#include "version.hpp"

////////////////// PRINT OUT CANDOCK CONFIGURATION FILE ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                help::Options::set_options(new Program::CmdLnOpts(argc, argv));

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info() <<
                        help::Options::get_options()->configuration_file() << endl;
                
        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
