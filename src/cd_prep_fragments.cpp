#include <iostream>
#include "candock/drm/drm.hpp"
#include "candock/program/cmdlnopts.hpp"
#include "candock/program/fragmentligands.hpp"
#include "version.hpp"

// FRAGMENTING OF LIGANDS

using namespace std;
using namespace candock;

int main(int argc, char* argv[]) {
    try {
        if (!drm::check_drm(Version::get_install_path() + "/.candock")) {
            throw logic_error(
                "CANDOCK has expired. Please contact your CANDOCK distributor "
                "to get a new version.");
        }

        help::Options::set_options(new Program::CmdLnOpts(
            argc, argv,
            Program::CmdLnOpts::STARTING | Program::CmdLnOpts::LIG_FRAMGENT));

        Benchmark main_timer;
        main_timer.display_time("Starting");

        cout << Version::get_banner() << Version::get_version()
             << Version::get_run_info();
        cout << help::Options::get_options()->configuration_file() << endl;

        Program::FragmentLigands ligand_fragmenter;
        ligand_fragmenter.run_step();

        main_timer.display_time("Finished");

    } catch (exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
    return 0;
}
