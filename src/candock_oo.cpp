#include <iostream>
#include <exception>
#include <typeinfo>

#include "program/cmdlnopts.hpp"
#include "program/findcentroids.hpp"
#include "program/fragmentligands.hpp"
#include "program/dockfragments.hpp"
#include "program/linkfragments.hpp"

#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecules.hpp"
#include "score/score.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "modeler/forcefield.hpp"
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

		//TODO: Combine into one class?????
		Program::Target targets (cmdl.get_string_option("target_dir"));
		targets.find_centroids(cmdl);
		
		Program::Target antitargets(cmdl.get_string_option("antitarget_dir"));
		antitargets.find_centroids(cmdl);

		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		//TODO: Make specific to each protein
		score.define_composition(targets.get_idatm_types(antitargets.get_idatm_types()),
								 ligand_fragmenter.ligand_idatm_types())
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		dbgmsg("START SCORE" << endl << score << "END SCORE");

		    targets.dock_fragments(score, ligand_fragmenter, cmdl);
		antitargets.dock_fragments(score, ligand_fragmenter, cmdl);

		cout << "Determining the best seeds to add" << endl;
		multiset<string>  target_seeds =     targets.determine_overlapping_seeds(cmdl.get_int_option("seeds_to_add"),   cmdl.get_int_option("seeds_till_good"));
		multiset<string> atarget_seeds = antitargets.determine_overlapping_seeds(cmdl.get_int_option("seeds_to_avoid"), cmdl.get_int_option("seeds_till_bad"));

		set<string> solo_target_seeds;
		std::set_difference( target_seeds.begin(),  target_seeds.end(),
							atarget_seeds.begin(), atarget_seeds.end(),
							std::inserter(solo_target_seeds, solo_target_seeds.end())
		);

		OMMIface::SystemTopology::loadPlugins();

		if (cmdl.get_bool_option("target_linking"))
			targets.link_fragments(score, cmdl);

		if (cmdl.get_bool_option("antitarget_linking"))
			antitargets.link_fragments(score,cmdl);

		Molib::PDBreader lpdb2(cmdl.prep_file(), Molib::PDBreader::all_models, 1);
		Molib::Molecules mol;
		lpdb2.parse_molecule(mol);

		cout << "Starting Design with " << solo_target_seeds.size() << " seeds." << endl;
		
		design::Design designer( mol.first() );
		designer.add_fragments_to_existing_molecule(common::read_top_seeds_files(solo_target_seeds, "targets/syk/" + cmdl.top_seeds_dir(), cmdl.top_seeds_file(), cmdl.top_percent() ));
		inout::output_file(designer.get_internal_designs(), "interal_designs.pdb");
		inout::output_file(designer.get_prepared_designs(), "designed.pdb");

		cmdl.display_time("finished");

	} catch ( exception& e) {
		cerr << e.what() << endl;
		return 1;
	}

	return 0;
}

