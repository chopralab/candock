#include <iostream>
#include "statchem/helper/logger.hpp"
#include "candock/program/cmdlnopts.hpp"
#include "version.hpp"

// PRINT OUT CANDOCK CONFIGURATION FILE

using namespace std;
using namespace candock;

int main(int argc, char* argv[]) {
    try {
        /*if (!drm::check_drm(Version::get_install_path() + "/.candock")) {
            throw logic_error(
                "CANDOCK has expired. Please contact your CANDOCK distributor "
                "to get a new version.");
        }*/

        statchem::Logger::set_all_stderr(true);

        help::Options::set_options(new Program::CmdLnOpts(argc, argv));

        cout << help::Options::get_options()->configuration_file() << endl;

    } catch (exception& e) {
        std::cerr << e.what() << endl;
        return 1;
    }
    return 0;
}
