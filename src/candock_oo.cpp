#include <iostream>
#include <exception>
#include <typeinfo>

#include "opts_candock.hpp"
#include "findcentroidsstep.hpp"
#include "fragmentligandsstep.hpp"
#include "dockfragmentsstep.hpp"
#include "linkfragmentsstep.hpp"

#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecules.hpp"
#include "score/score.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "modeler/forcefield.hpp"
#include "modeler/systemtopology.hpp"

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
 * |->---------------------> Dock Fragments <------------------------|
 * |                               |                                 |
 * |                               V                                 |
 * L-> --------------------> Link Framgents <-------------------------
 * 
 * **************************************************************************/

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;

		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

		/* Run section of Candock designed to find binding sites
 		 * Currently, this runs ProBIS and does not require any
 		 * previous step to be competed.
 		 *
		 */

		Program::FindCentroidsStep find_centroids(receptors[0]);
		find_centroids.run_step(cmdl);

		Program::FragmentLigandsStep ligand_fragmenter;
		ligand_fragmenter.run_step(cmdl);

		/* Create receptor grid
		 * 
		 */
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());

		/* Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_forcefield_file(cmdl.amber_xml_file())
			.parse_forcefield_file(cmdl.water_xml_file());

		receptors[0].prepare_for_mm(ffield, gridrec);

		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		score.define_composition(receptors.get_idatm_types(), ligand_fragmenter.ligand_idatm_types())
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		dbgmsg("START SCORE" << endl << score << "END SCORE");

		Program::DockFragmentsStep fragment_docker( find_centroids, ligand_fragmenter,
													score, gridrec );
		
		/**
		 * Insert topology for cofactors, but not for standard residues
		 * that are already known to forcefield (e.g., amino acid residues)
		 *
		 */
		ffield.insert_topology(receptors[0]); // fixes issue #115

		OMMIface::SystemTopology::loadPlugins();

		Program::LinkFragmentsStep link_fragments( receptors[0], score, ffield, gridrec);
		link_fragments.run_step(cmdl);
		
	} catch ( exception& e) {
		cerr << e.what() << endl;
	}
}
