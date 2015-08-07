#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "opts_candock.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "common.hpp"
#include "helper/error.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/pdbreader.hpp"
#include "openmm/forcefield.hpp"
#include "openmm/moleculeinfo.hpp"
#include "score/score.hpp"
#include "openmm/omm.hpp"
#include "linker/linker.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"
using namespace std;

CmdLnOpts cmdl;

////////////////// TEST MINIMIZATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		inout::output_file("", cmdl.mini_ligands_file()); // output docked & minimized ligands conformations

		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model|Molib::PDBreader::skip_hetatm|Molib::PDBreader::hydrogens);
		Molib::Molecules receptors = rpdb.parse_molecule();

		// if empty, add dummy receptor
		if (receptors.empty())
			receptors.add(new Molib::Molecule("dummy"));

		receptors[0].filter(Molib::Residue::protein, cmdl.receptor_chain_id());

			
		Molib::PDBreader lpdb(cmdl.ligand_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

		/* Compute atom types for receptor (gaff types not needed since 
		 * they are read from the forcefield xml file)
		 * 
		 */
		receptors.compute_idatm_type();

		/* Create receptor grid
		 * 
		 */
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());
		

		Molib::Molecules ligands = lpdb.parse_molecule();

		// if empty, add dummy receptor
		if (ligands.empty())
			ligands.add(new Molib::Molecule("dummy"));

		set<int> ligand_idatm_types;
		ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);


		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			gridrec, cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond());

		/** 
		 * Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond(), cmdl.scale_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file());

	
		receptors[0].prepare_for_mm(ffield, gridrec);

		OMMIface::SystemTopology::loadPlugins();

		for (auto &ligand : ligands) {
			ffield.insert_topology(ligand);
			try {
				/**
				 * Minimize system
				 */
				 
				OMMIface::Modeler modeler;
				modeler.set_forcefield(ffield);
				modeler.set_forcefield_type(cmdl.fftype());
				modeler.set_distance_cutoff(cmdl.dist_cutoff());
				modeler.set_use_constraints(false);
				modeler.set_step_size_in_fs(2.0);
				modeler.set_tolerance(cmdl.tolerance());
				modeler.set_max_iterations(cmdl.max_iterations());
				modeler.set_update_freq(cmdl.update_freq());
				modeler.set_position_tolerance(cmdl.position_tolerance());
				
				modeler.add_topology(receptors[0].get_atoms());
				modeler.add_topology(ligand.get_atoms());
				
				modeler.add_crds(receptors[0].get_atoms(), receptors[0].get_crds());
				modeler.add_crds(ligand.get_atoms(), ligand.get_crds());
				
				modeler.init_openmm();
				
				modeler.unmask(receptors[0].get_atoms());
				modeler.unmask(ligand.get_atoms());
				
				modeler.minimize_state();
				
				Molib::Molecule mini_receptor = modeler.get_state(receptors[0].get_atoms());
				Molib::Molecule mini_ligand = modeler.get_state(ligand.get_atoms());
				
				mini_state.undo_mm_specific();
				inout::output_file(mini_state, cmdl.mini_ligands_file(), ios_base::app);

				
				//~ OMMIface::OMM omm(receptors[0], ligand, ffield, 
					//~ cmdl.fftype(), cmdl.dist_cutoff());
				//~ omm.minimize(cmdl.tolerance(), cmdl.max_iterations(), cmdl.update_freq()); // minimize
				//~ dbgmsg(ligand);
				//~ auto ret = omm.get_state(receptors[0], ligand);
				//~ Molib::Molecules mini;
				//~ Molib::Molecule &minimized_receptor = 
					//~ mini.add(new Molib::Molecule(ret.first));
				//~ Molib::Molecule &minimized_ligand = 
					//~ mini.add(new Molib::Molecule(ret.second));
		//~ 
				//~ minimized_receptor.undo_mm_specific();
		//~ 
				//~ dbgmsg("MOLECULES AFTER MINIMIZATION : " << endl 
					//~ << minimized_ligand << endl 
					//~ << minimized_receptor);
				//~ inout::output_file(mini, cmdl.mini_ligands_file(), ios_base::app);

				//~ OMMIface::Energies energies;
				//~ energies[&minimized_ligand] = omm.get_energy_components(
					//~ minimized_receptor, minimized_ligand, cmdl.dist_cutoff());
				//~ minimized_receptor.undo_mm_specific();
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
