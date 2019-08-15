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
#include <exception>
#include <typeinfo>

#include <boost/filesystem.hpp>

#include "program/cmdlnopts.hpp"
#include "program/target.hpp"
#include "program/fragmentligands.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

using namespace std;

/*****************************************************************************
 *
 * <---receptor.pdb                                             ligands.mol2
 * |        |                                                        |
 * |        |                                                        |
 * |        V                                                        V
 * |   Find Centroids                                         Fragment Ligands
 * |        |                                                        |
 * |        V                                                        |
 * |->------>--------------> Dock Fragments <----------------------<-|
 * |                               |                                 |
 * |                               V                                 |
 * L-> --------------------> Link Framgents <----------------------<-|
 *                                 |                                 |
 *                                 V                                 |
 *                        Design of Compounds                        |
 * 
 * **************************************************************************/

int main(int argc, char* argv[]) {
        try {

                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                help::Options::set_options( new Program::CmdLnOpts( 
                    argc, argv, Program::CmdLnOpts::DESIGN));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                cout << help::Options::get_options()->configuration_file() << endl;

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                //TODO: Combine into one class?????
                Program::Target targets (cmdl.get_string_option("target_dir"));
                targets.link_fragments(ligand_fragmenter);

                Program::Target antitargets(cmdl.get_string_option("antitarget_dir"));
                antitargets.find_centroids();
                antitargets.dock_fragments(ligand_fragmenter);

                if (cmdl.get_bool_option("antitarget_linking"))
                        antitargets.link_fragments(ligand_fragmenter);

                set<string> solo_target_seeds = Program::Target::determine_non_overlapping_seeds(targets, antitargets);
                if (cmdl.get_bool_option("new_scaffold") || ! boost::filesystem::is_regular_file(cmdl.get_string_option("prep"))) {
                        targets.make_scaffolds(ligand_fragmenter, solo_target_seeds);
                }

                targets.design_ligands(ligand_fragmenter, solo_target_seeds);

                main_timer.display_time("Finished");
        } catch ( exception& e) {
                cerr << e.what() << endl;
                return 1;
        }

        return 0;
}

