/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include <iostream>
#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "version.hpp"
#include "drm/drm.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {
                if( !drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                help::Options::set_options(new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::LIG_FRAMGENT ));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
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
