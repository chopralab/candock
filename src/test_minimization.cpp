#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include <boost/filesystem.hpp>
#include "program/cmdlnopts.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "helper/path.hpp"
#include "program/common.hpp"
#include "helper/error.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/nrset.hpp"
#include "pdbreader/pdbreader.hpp"
#include "modeler/forcefield.hpp"
#include "score/score.hpp"
#include "modeler/modeler.hpp"
#include "linker/linker.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"

using namespace std;
using namespace Program;

////////////////// TEST MINIMIZATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv, CmdLnOpts::STARTING | CmdLnOpts::FORCE_FIELD | CmdLnOpts::SCORING);
		cmdl.display_time("started");
		cout << cmdl << endl;

		//inout::output_file("", cmdl.docked_file()); // output docked & minimized ligands conformations

		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

		// if empty, add dummy receptor
		if (receptors.empty())
			receptors.add(new Molib::Molecule("dummy"));

		Molib::PDBreader lpdb(cmdl.prep_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

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
		

		Molib::Molecules ligands = lpdb.parse_molecule();

		// if empty, add dummy receptor
		if (ligands.empty())
			ligands.add(new Molib::Molecule("dummy"));

		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		score.define_composition(receptors.get_idatm_types(), ligands.get_idatm_types())
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		dbgmsg("START SCORE" << endl << score << "END SCORE");

		/** 
		 * Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file())
			.parse_forcefield_file(cmdl.water_xml_file());

	
		receptors[0].prepare_for_mm(ffield, gridrec);

		/**
		 * Insert topology for cofactors, but not for standard residues
		 * that are already known to forcefield (e.g., amino acid residues)
		 *
		 */
		ffield.insert_topology(receptors[0]);

		OMMIface::SystemTopology::loadPlugins();

		for (auto &ligand : ligands) {
			try {
				
				/**
				 * Minimize system
				 */
				
				ffield.insert_topology(ligand);
				//~ ligand.set_name("org");

				OMMIface::Modeler modeler(ffield, cmdl.fftype(), cmdl.dist_cutoff(),
					cmdl.tolerance(), cmdl.max_iterations(), cmdl.update_freq(), cmdl.position_tolerance(),
					false, 2.0);
						
				modeler.add_topology(receptors[0].get_atoms());
				modeler.add_topology(ligand.get_atoms());

				modeler.init_openmm();
//~ 
				//~ // change coordinate of some ligand atoms
				//~ Molib::Atom::Vec substruct;
				//~ for (auto &patom : ligand.get_atoms()) {
					//~ if (patom->atom_number() == 2 || patom->atom_number() == 1 
						//~ || patom->atom_number() == 4 || patom->atom_number() == 57) {
						//~ substruct.push_back(patom);
					//~ }
				//~ }
				//~ for (auto &patom : substruct) {
					//~ patom->set_crd(patom->crd() + 5.0);
				//~ }

				modeler.add_crds(receptors[0].get_atoms(), receptors[0].get_crds());
				modeler.add_crds(ligand.get_atoms(), ligand.get_crds());
				
				modeler.init_openmm_positions();
				
				//~ modeler.unmask(receptors[0].get_atoms());
				//~ modeler.mask(ligand.get_atoms());
				//~ modeler.unmask(substruct);

				cout << "Initial energy = " << score.non_bonded_energy(gridrec, ligand) << endl;

				modeler.minimize_state(ligand, receptors[0], score);

				// init with minimized coordinates
				Molib::Molecule minimized_receptor(receptors[0], modeler.get_state(receptors[0].get_atoms()));
				Molib::Molecule minimized_ligand(ligand, modeler.get_state(ligand.get_atoms()));

				minimized_receptor.undo_mm_specific();

				Molib::Atom::Grid gridrec(minimized_receptor.get_atoms());
				const double energy = score.non_bonded_energy(gridrec, minimized_ligand);
				
				cout << "Minimized energy = " << energy << endl;

				inout::output_file(Molib::Molecule::print_complex(minimized_ligand, minimized_receptor, energy), 
					Path::join(cmdl.docked_dir(), minimized_ligand.name() + ".pdb" ), ios_base::app);
				
				//~ /**
				 //~ * This section is a test to see if you can change positions 
				//~ * of some atoms without reinitializing the whole openmm
				 //~ */
				 //~ 
				//~ // change coordinate of ligand atoms
				//~ ligand.set_name("mod");
				//~ for (auto &patom : ligand.get_atoms()) {
					//~ patom->set_crd(patom->crd() + 2.0);
				//~ }
				//~ // initialize only new ligand coordinates
				//~ modeler.add_crds(ligand.get_atoms(), ligand.get_crds());
				//~ modeler.init_openmm_positions();
				//~ 
				//~ modeler.minimize_state(ligand, receptors[0], score);
//~ 
				//~ // init with minimized coordinates
				//~ Molib::Molecule mod_minimized_receptor(receptors[0], modeler.get_state(receptors[0].get_atoms()));
				//~ Molib::Molecule mod_minimized_ligand(ligand, modeler.get_state(ligand.get_atoms()));
//~ 
				//~ mod_minimized_receptor.undo_mm_specific();
//~ 
				//~ Molib::Atom::Grid mod_gridrec(mod_minimized_receptor.get_atoms());
				//~ const double mod_energy = score.non_bonded_energy(mod_gridrec, mod_minimized_ligand);
//~ 
				//~ inout::output_file(Molib::Molecule::print_complex(mod_minimized_ligand, mod_minimized_receptor, mod_energy), 
					//~ "mod_" + cmdl.mini_file(), ios_base::app);
					//~ 
				//~ /*****************************/

			} catch (exception& e) {
				cerr << "MINIMIZATION FAILED FOR LIGAND " << ligand.name() 
					<< " because of " << e.what() << endl;
			}
			ffield.erase_topology(ligand);

		}
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
