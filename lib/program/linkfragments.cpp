#include "linkfragments.hpp"

#include "common.hpp"
#include "linker/linker.hpp"
#include "modeler/modeler.hpp"
#include "helper/path.hpp"
#include "helper/grep.hpp"

namespace Program {

	bool LinkFragments::__can_read_from_files ( const CmdLnOpts& cmdl )
	{
		boost::regex regex;
		regex.assign("REMARK   5 MOLECULE (\\w*)");
		std::ifstream file( cmdl.prep_file() );
		
		std::vector<std::string> all_names = Grep::search_stream(file, regex);
		
		boost::filesystem::path p ( __receptor.name());
		p /= cmdl.docked_dir();
		
		bool is_done = true;
		
		for (auto molec : all_names) {
			boost::filesystem::path p2 = p / (molec + ".pdb");
			if ( inout::Inout::file_size(p2.string()) <= 0) {
				is_done = false;
			} else {
				Molib::PDBreader conf(p2.string(), Molib::PDBreader::skip_atom | Molib::PDBreader::first_model, 1);
				conf.parse_molecule(__all_top_poses);
			}
		}
		
		return is_done;
	}

	void LinkFragments::__read_from_files ( const CmdLnOpts& cmdl )
	{
		cout << "Linking for all molecules in " << cmdl.prep_file() << " for " << __receptor.name() << " is complete, skipping." << endl;
	}

	void LinkFragments::__link_ligand( Molib::Molecule ligand, const CmdLnOpts& cmdl, const OMMIface::ForceField& ffield ) {
		boost::filesystem::path p(__receptor.name());
		p = p / cmdl.docked_dir() / (ligand.name() + ".pdb");
		if ( inout::Inout::file_size(p.string()) > 0) {
			cout << ligand.name() << " is alread docked to " << __receptor.name() << ", reading from already docked file." << endl;
			return;
		}

		// if docking of one ligand fails, docking of others shall continue...
		try { 
			dbgmsg("LINKING LIGAND : " << endl << ligand);

			if ( ligand.first().first().get_rigid().size() < 1 ) {
				throw Error("No seeds to link");
			}

			/** 
				* Read top seeds for this ligand
				*/
			Molib::NRset top_seeds = common::read_top_seeds_files(ligand,
					Path::join(__receptor.name(), cmdl.top_seeds_dir()), cmdl.top_seeds_file(), cmdl.top_percent());

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
			OMMIface::Modeler modeler(ffield, cmdl.fftype(), cmdl.dist_cutoff(),
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
			Linker::DockedConformation::sort(docks);

			int model = 0;
			for (auto &docked : docks) {
				common::change_residue_name(docked.get_ligand(), "CAN"); 
				inout::output_file(Molib::Molecule::print_complex(docked.get_ligand(), docked.get_receptor(), docked.get_energy(), ++model), 
						p.string(), ios_base::app); // output docked molecule conformations
			}
			__all_top_poses.add( new Molib::Molecule (docks[0].get_ligand()) );
		} catch (exception& e) {
			cerr << "Error: skipping ligand " << ligand.name() << " with " << __receptor.name() << " due to : " << e.what() << endl;
			stringstream ss;
			common::change_residue_name(ligand, "CAN");
			ss << "REMARK  20 non-binder " << ligand.name() << " with " << __receptor.name() << " because " << e.what() << endl << ligand;
			inout::Inout::file_open_put_stream(p.string(), ss, ios_base::app);
		} 
	}

	void LinkFragments::__continue_from_prev ( const CmdLnOpts& cmdl )
	{

		cout << "Starting to dock the fragments into originally given ligands" << endl;

		Molib::PDBreader lpdb2(cmdl.prep_file(), Molib::PDBreader::all_models, 1);

		std::vector<std::thread> threads;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back( std::thread([&,this] {
				OMMIface::ForceField ffcopy(__ffield);
				Molib::Molecules ligands;

				while (lpdb2.parse_molecule(ligands)) {
					Molib::Molecule &ligand = ligands.first();

					/**
					 * Ligand's resn MUST BE UNIQUE for ffield
					 */
					common::change_residue_name(ligand, __concurrent_numbering, __ligand_cnt);
					ffcopy.insert_topology(ligand);
					__link_ligand(ligand, cmdl, ffcopy);
					ffcopy.erase_topology(ligand); // he he
					ligands.clear();
				}
			} ) );
		}

		for(auto& thread : threads) {
			thread.join();
		}

		cout << "Linking of fragments is complete" << endl;
	}

	void LinkFragments::link_ligands(const Molib::Molecules ligands, const CmdLnOpts& cmdl) {
		int j = 0;

		std::vector<std::thread> threads;
		std::mutex counter_protect;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back( std::thread([&,this] {
				OMMIface::ForceField ffcopy(__ffield);

				while (true)
				{
					unique_lock<std::mutex> guard(counter_protect, std::defer_lock);
					guard.lock();

					if( j >= ligands.size() )
						return;

					Molib::Molecule &ligand = ligands[j++];
					
					guard.unlock();

					/**
					 * Ligand's resn MUST BE UNIQUE for ffield
					 */
					common::change_residue_name(ligand, __concurrent_numbering, __ligand_cnt);
					ffcopy.insert_topology(ligand);
					__link_ligand(ligand, cmdl, ffcopy);
					ffcopy.erase_topology(ligand); // he he
				}
			} ) );
		}

		for(auto& thread : threads) {
			thread.join();
		}
	}
}
