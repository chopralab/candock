#include <iostream>
#include <exception>
#include <typeinfo>

#include "opts_candock.hpp"
#include "fragmentligands.hpp"

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

		Program::Target targets (cmdl.get_string_option("target_dir"));
		targets.find_centroids(cmdl);
		
		//Program::Target anitargets(cmdl, "antitarget_dir");

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

		targets.determine_overlapping_seeds(cmdl.get_int_option("seeds_to_add"), cmdl.get_int_option("seeds_till_good"), false);

		OMMIface::SystemTopology::loadPlugins();
		targets.link_fragments(score, cmdl);

		cmdl.display_time("finished");

	} catch ( exception& e) {
		cerr << e.what() << endl;
		return 1;
	}

	return 0;
}

