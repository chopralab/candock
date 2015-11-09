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
#include "score/score.hpp"
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
		
		inout::output_file("", cmdl.docked_file()); // output docked & minimized ligands conformations

		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model|Molib::PDBreader::skip_hetatm|Molib::PDBreader::hydrogens);
		Molib::Molecules receptors = rpdb.parse_molecule();

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

		set<int> ligand_idatm_types;
		ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);


		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond(),
			cmdl.scale_non_bond());
			
		dbgmsg(score);

		for (auto &ligand : ligands) {

			const double energy = score.non_bonded_energy(gridrec, ligand);
			
			inout::output_file(Molib::Molecule::print_complex(ligand, receptors[0], energy), 
				cmdl.docked_file(), ios_base::app);
		}
		
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
