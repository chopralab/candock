#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "program/cmdlnopts.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "helper/path.hpp"
#include "helper/error.hpp"
#include "molib/grid.hpp"
#include "molib/molecules.hpp"
#include "molib/pdbreader.hpp"
#include "score/score.hpp"

using namespace std;
using namespace Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv, CmdLnOpts::STARTING | CmdLnOpts::SCORING);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		//inout::output_file("", cmdl.docked_file()); // output docked & minimized ligands conformations

		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

		/** 
		 * Compute atom types for receptor and cofactor(s): gaff types for protein, 
		 * Mg ions, and water are read from the forcefield xml file later on while 
		 * gaff types for cofactors (ADP, POB, etc.) are calculated de-novo here
		 * 
		 */
		receptors.compute_idatm_type()
			.compute_hydrogen()
			.compute_bond_order()
			.compute_bond_gaff_type()
			.refine_idatm_type()
			.erase_hydrogen()  // needed because refine changes connectivities
			.compute_hydrogen()   // needed because refine changes connectivities
			.compute_ring_type()
			.compute_gaff_type()
			.compute_rotatable_bonds() // relies on hydrogens being assigned
			.erase_hydrogen();

		/* Create receptor grid
		 * 
		 */
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());
		

		Molib::PDBreader lpdb(cmdl.prep_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

		Molib::Molecules ligands = lpdb.parse_molecule();

		set<int> ligand_idatm_types;
		ligand_idatm_types = ligands.get_idatm_types(ligand_idatm_types);

		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		score.define_composition(receptors.get_idatm_types(), ligand_idatm_types)
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		for (auto &ligand : ligands) {

			const double energy = score.non_bonded_energy(gridrec, ligand);
			
			inout::output_file(Molib::Molecule::print_complex(ligand, receptors[0], energy), 
				Path::join( cmdl.docked_dir(), ligand.name() + ".pdb"),  ios_base::app);
		}
		
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
