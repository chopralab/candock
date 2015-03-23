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

////////////////// TEST FRAGMENTING OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model|Molib::PDBreader::skip_hetatm);
		//~ Molib::Molecules receptors = rpdb.parse_molecule();
		
		Molib::PDBreader lpdb(cmdl.ligand_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

		//~ Molib::Molecules ligands = lpdb.parse_molecule();

		Molib::Molecules seeds;
		set<string> added;
		while(1 != 0) {
			Molib::Molecules ligands = lpdb.parse_molecule();
			if (ligands.empty()) break;

			// Compute properties, such as idatm atom types, fragments, seeds,
			// rotatable bonds
			//~ ligands.compute_idatm_type()
				//~ .compute_ring_type()
				//~ .compute_rotatable_bonds()
				//~ .compute_overlapping_rigid_segments()
				//~ .compute_seeds(cmdl.seeds_file());
			ligands.compute_idatm_type()
				.compute_hydrogen()
				.compute_bond_order()
				.compute_ring_type()
				.compute_bond_gaff_type()
				.compute_gaff_type()
				.compute_rotatable_bonds() // relies on hydrogens being assigned
				.erase_hydrogen()
				.compute_overlapping_rigid_segments()
				.compute_seeds(cmdl.seeds_file());
			//~ common::create_mols_from_seeds(added, seeds, ligands);
			common::create_mols_from_fragments(added, seeds, ligands);
			//~ inout::output_file(ligands, cmdl.prep_file(), ios_base::app);
		}
		int i = 0;
		for (auto &seed : seeds) {
			inout::output_file(seed, cmdl.prep_file() + "." + help::to_string(i++));
		}
		dbgmsg("SEEDS FOUND IN THE MOLECULE : " << endl << seeds);
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
