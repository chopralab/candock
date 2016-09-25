#include <iostream>
#include <exception>
#include <typeinfo>

#include "opts_candock.hpp"
#include "findcentroids.hpp"
#include "fragmentligands.hpp"
#include "dockfragments.hpp"
#include "linkfragments.hpp"

#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecules.hpp"
#include "score/score.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "modeler/forcefield.hpp"
#include "modeler/systemtopology.hpp"
#include "helper/inout.hpp"
#include "target.hpp"

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
 * L-> --------------------> Link Framgents <-------------------------
 * 
 * **************************************************************************/

int main(int argc, char* argv[]) {
	try {
		Program::CmdLnOpts cmdl;
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		Program::FragmentLigands ligand_fragmenter;
		ligand_fragmenter.run_step(cmdl);

		Program::Target targets(cmdl, "target_dir");

		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		score.define_composition(targets.get_idatm_types(), ligand_fragmenter.ligand_idatm_types())
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		dbgmsg("START SCORE" << endl << score << "END SCORE");

		targets.dock_fragments(score, ligand_fragmenter, cmdl);

// 		OMMIface::ForceField ffield;
// 		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
// 			.add_kb_forcefield(score, cmdl.step_non_bond())
// 			.parse_forcefield_file(cmdl.amber_xml_file())
// 			.parse_forcefield_file(cmdl.water_xml_file());
// 
// 		receptors[0].prepare_for_mm(ffield, gridrec);
// 		
// 		/**
// 		 * Insert topology for cofactors, but not for standard residues
// 		 * that are already known to forcefield (e.g., amino acid residues)
// 		 *
// 		 */
// 		ffield.insert_topology(receptors[0]); // fixes issue #115
// 
// 		OMMIface::SystemTopology::loadPlugins();
// 
// 		Program::LinkFragments link_fragments( receptors[0], score, ffield, gridrec);
// 		link_fragments.run_step(cmdl);
// 
// 		cmdl.display_time("finished");

	} catch ( exception& e) {
		cerr << e.what() << endl;
		return 1;
	}

	return 0;
}

