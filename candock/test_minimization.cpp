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
//~ #include "pdbreader/idatm.hpp"
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
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model|Molib::PDBreader::skip_hetatm);
		Molib::Molecules receptors = rpdb.parse_molecule();
		
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
		Molib::MolGrid gridrec(receptors[0].get_atoms(cmdl.receptor_chain_id(), 
			Molib::Residue::protein));

		//~ {
			//~ OMMIface::ForceField ffield;
			//~ ffield
				//~ .parse_forcefield_file(cmdl.amber_xml_file());
			//~ receptors[0].prepare_for_mm(ffield, gridrec);
			//~ throw Error ("exit after prepare for mm");
		//~ }

		Molib::Molecules ligands = lpdb.parse_molecule();

		ligands.compute_idatm_type();
		//~ throw Error("exit after idatm typing");
		ligands.compute_hydrogen();
		//~ throw Error("exit after compute hydrogen");
		ligands.compute_bond_order();
		//~ throw Error("exit after compute bond order");
		ligands.compute_ring_type();
		//~ throw Error("exit after ring typing");
		ligands.compute_bond_gaff_type();
		//~ throw Error("exit after bond gaff typing");
		ligands.compute_gaff_type();
		//~ throw Error("exit after ligand atom typing");
		ligands.erase_hydrogen();
		//~ throw Error("exit after erasing hydrogens");
		ligands.compute_rotatable_bonds();
		//~ throw Error("exit after computing rotatable bonds");
		ligands.compute_overlapping_rigid_segments();
		//~ throw Error("exit after computing overlapping rigid segments");
		
		set<int> ligand_idatm_types;
		ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);


		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			gridrec, cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond());

		/* Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_residue_topology(ligands)
			.add_kb_forcefield(score, cmdl.step_non_bond(), cmdl.scale_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file());

		receptors[0].prepare_for_mm(ffield, gridrec);
			
		OMMIface::OMM omm(receptors[0], ligands.first(), ffield, 
			cmdl.fftype(), cmdl.dist_cutoff());
		omm.minimize(cmdl.tolerance(), cmdl.max_iterations()); // minimize
		//~ omm.minimize(0.0000000000000001, cmdl.max_iterations()); // minimize
		dbgmsg(ligands.first());
		auto ret = omm.get_state(receptors[0], ligands.first());
		Molib::Molecules mini;
		Molib::Molecule &minimized_receptor = 
			mini.add(new Molib::Molecule(ret.first));
		Molib::Molecule &minimized_ligand = 
			mini.add(new Molib::Molecule(ret.second));

		minimized_receptor.undo_mm_specific();

		dbgmsg("MOLECULES AFTER MINIMIZATION : " << endl 
			<< minimized_ligand << endl 
			<< minimized_receptor);
		//~ OMMIface::Energies energies;
		//~ energies[&minimized_ligand] = omm.get_energy_components(
			//~ minimized_receptor, minimized_ligand, cmdl.dist_cutoff());
		//~ minimized_receptor.undo_mm_specific();
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
