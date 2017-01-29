#include <iostream>
#include <exception>
#include <typeinfo>

#include "program/cmdlnopts.hpp"
#include "program/findcentroids.hpp"
#include "program/fragmentligands.hpp"
#include "program/dockfragments.hpp"
#include "program/linkfragments.hpp"

#include "pdbreader/molecules.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "modeler/systemtopology.hpp"
#include "helper/inout.hpp"
#include "program/target.hpp"
#include "design/design.hpp"

#include <algorithm>
#include "program/common.hpp"

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
                Program::CmdLnOpts *cmdlopts = new Program::CmdLnOpts;
                cmdlopts->init(argc, argv, Program::CmdLnOpts::DESIGN);

                Benchmark main_timer;
                main_timer.display_time("started");
                cout << *cmdlopts << endl;

                help::Options::set_options(cmdlopts);

                Program::FragmentLigands ligand_fragmenter;
                ligand_fragmenter.run_step();

                //TODO: Combine into one class?????
                Program::Target targets (cmdl.get_string_option("target_dir"));
                targets.find_centroids();
                targets.dock_fragments(ligand_fragmenter);

                Program::Target antitargets(cmdl.get_string_option("antitarget_dir"));
                antitargets.find_centroids();
                antitargets.dock_fragments(ligand_fragmenter);

                OMMIface::SystemTopology::loadPlugins();

                targets.link_fragments();

                if (cmdl.get_bool_option("antitarget_linking"))
                        antitargets.link_fragments();

                set<string> solo_target_seeds = Program::Target::determine_non_overlapping_seeds(targets, antitargets);

                targets.design_ligands(ligand_fragmenter, solo_target_seeds);

                main_timer.display_time("finished");
        } catch ( exception& e) {
                cerr << e.what() << endl;
                return 1;
        }

        return 0;
}

