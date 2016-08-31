#include "linkfragmentsstep.hpp"

#include "common.hpp"
#include "linker/linker.hpp"
#include "modeler/modeler.hpp"
#include "helper/path.hpp"

namespace Program {

	bool LinkFragmentsStep::__can_read_from_files ( const CmdLnOpts& cmdl )
	{
		return false;
	}

	void LinkFragmentsStep::__read_from_files ( const CmdLnOpts& cmdl )
	{
		
	}

	void Program::LinkFragmentsStep::__link_ligand_from_fragment( Molib::PDBreader& lpdb2, const CmdLnOpts& cmdl)
	{
		OMMIface::ForceField ffcopy(__ffield);
		Molib::Molecules ligands;

		while (lpdb2.parse_molecule(ligands)) {
	
			Molib::Molecule &ligand = ligands.first();
			dbgmsg("LINKING LIGAND : " << endl << ligand);

			/**
			 * Ligand's resn MUST BE UNIQUE for ffield
			 */
			common::change_residue_name(ligand, __concurrent_numbering, __ligand_cnt); 

			// if docking of one ligand fails, docking of others shall continue...
			try { 

				ffcopy.insert_topology(ligand);

				/** 
				 * Read top seeds for this ligand
				 */	 
				Molib::NRset top_seeds = common::read_top_seeds_files(ligand,
						cmdl.top_seeds_dir(), cmdl.top_seeds_file(), cmdl.top_percent());
	
				ligand.erase_properties(); // required for graph matching
				top_seeds.erase_properties(); // required for graph matching
	
				/** 
				 * Jiggle the coordinates by one-thousand'th of an Angstrom to avoid minimization failures
				 * with initial bonded relaxation failed errors
				 */
				top_seeds.jiggle();
	
				/* Init minization options and constants, including ligand and receptor topology
				 *
				 */
				OMMIface::Modeler modeler(ffcopy, cmdl.fftype(), cmdl.dist_cutoff(),
						cmdl.tolerance(), cmdl.max_iterations(), cmdl.update_freq(), 
						cmdl.position_tolerance(), false, 2.0);

				/**
				 * Connect seeds with rotatable linkers, account for symmetry, optimize 
				 * seeds with appendices, minimize partial conformations between linking.
				 * 
				 */
				Linker::Linker linker(modeler, __receptor, ligand, top_seeds, __gridrec, __score, 
						cmdl.iterative(), cmdl.dist_cutoff(), cmdl.spin_degrees(), 
						cmdl.tol_seed_dist(), cmdl.lower_tol_seed_dist(), 
						cmdl.upper_tol_seed_dist(), cmdl.max_possible_conf(),
						cmdl.link_iter(), cmdl.clash_coeff(), cmdl.docked_clus_rad(), 
						cmdl.max_allow_energy(), cmdl.max_num_possibles(),
						cmdl.max_clique_size(), cmdl.max_iterations_final());

				Linker::DockedConformation::Vec docks = linker.link();

				for (auto &docked : docks) {
					common::change_residue_name(docked.get_ligand(), "CAN"); 
					inout::output_file(Molib::Molecule::print_complex(docked.get_ligand(), docked.get_receptor(), docked.get_energy()), 
							Path::join(cmdl.docked_dir(), ligand.name() + ".pdb" ) , ios_base::app); // output docked molecule conformations
				}
			} catch (exception& e) { 
				cerr << "Error: skipping ligand " << ligand.name() << " due to : " << e.what() << endl; 
			} 
			ffcopy.erase_topology(ligand); // he he
			ligands.clear();
		}
	}

	void LinkFragmentsStep::__continue_from_prev ( const CmdLnOpts& cmdl )
	{
		Molib::PDBreader lpdb2(cmdl.ligand_file(), Molib::PDBreader::all_models, 1);

		std::vector<std::thread> threads;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back( std::thread([&,this] {__link_ligand_from_fragment(lpdb2, cmdl);} ) );
		}
		for(auto& thread : threads) {
			thread.join();
		}
	}

}
