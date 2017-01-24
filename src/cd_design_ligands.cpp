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
		Program::CmdLnOpts cmdl;
		cmdl.init(argc, argv, Program::CmdLnOpts::DESIGN);
		cmdl.display_time("started");
		cout << cmdl << endl;

		Program::FragmentLigands ligand_fragmenter;
		ligand_fragmenter.run_step(cmdl);

		//TODO: Combine into one class?????
		Program::Target targets (cmdl.get_string_option("target_dir"));
		targets.find_centroids(cmdl);
		targets.dock_fragments(ligand_fragmenter, cmdl);

		Program::Target antitargets(cmdl.get_string_option("antitarget_dir"));
		antitargets.find_centroids(cmdl);
		antitargets.dock_fragments(ligand_fragmenter, cmdl);

		OMMIface::SystemTopology::loadPlugins();
                
		targets.link_fragments(cmdl);

		if (cmdl.get_bool_option("antitarget_linking"))
			antitargets.link_fragments(cmdl);

                
                set<string> solo_target_seeds = Program::Target::determine_non_overlapping_seeds(targets, antitargets, cmdl);
                cout << "ASDFASDF " <<  inout::Inout::file_size(cmdl.prep_file()) << endl;
                if (/*cmdl.get_bool_option("new_scaffold") || */ ! inout::Inout::file_size(cmdl.prep_file())) {
                    cout << "Good" << endl;
                        targets.make_scaffolds(cmdl, ligand_fragmenter, solo_target_seeds);
                }

		targets.design_ligands(cmdl, ligand_fragmenter, solo_target_seeds);

		cmdl.display_time("finished");
	} catch ( exception& e) {
		cerr << e.what() << endl;
		return 1;
	}

	return 0;
}

