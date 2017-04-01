#include <iostream>
#include <exception>
#include <typeinfo>

#include "program/cmdlnopts.hpp"
#include "program/target.hpp"
#include "program/fragmentligands.hpp"

#include "version.hpp"

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
                help::Options::set_options( new Program::CmdLnOpts( 
                    argc, argv, Program::CmdLnOpts::DESIGN));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info() <<
                        help::Options::get_options()->configuration_file() << endl;

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                //TODO: Combine into one class?????
                Program::Target targets (cmdl.get_string_option("target_dir"));
                targets.find_centroids();
                targets.dock_fragments(ligand_fragmenter);

                Program::Target antitargets(cmdl.get_string_option("antitarget_dir"));
                antitargets.find_centroids();
                antitargets.dock_fragments(ligand_fragmenter);

                targets.link_fragments();

                if (cmdl.get_bool_option("antitarget_linking"))
                        antitargets.link_fragments();

                set<string> solo_target_seeds = Program::Target::determine_non_overlapping_seeds(targets, antitargets);
                if (/*cmdl.get_bool_option("new_scaffold") || */ ! boost::filesystem::is_regular_file(cmdl.get_string_option("prep"))) {
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

