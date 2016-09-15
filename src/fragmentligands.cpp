#include "fragmentligands.hpp"

#include <boost/filesystem.hpp>
#include <mutex>

#include "helper/inout.hpp"
#include "common.hpp"
#include "opts_candock.hpp"
#include "pdbreader/molecules.hpp"

namespace Program {
	
	bool FragmentLigands::__can_read_from_files( const CmdLnOpts& cmdl ) {
		return inout::Inout::file_size( cmdl.prep_file() ) > 0 && inout::Inout::file_size( cmdl.seeds_file() ) > 0;
	}

	void FragmentLigands::__read_from_files( const CmdLnOpts& cmdl) {

		cout << "Reading fragmented files from " << cmdl.prep_file() << endl;

		Molib::PDBreader lpdb(cmdl.prep_file(), Molib::PDBreader::all_models, cmdl.max_num_ligands());

		Molib::Molecules ligands;
		while(lpdb.parse_molecule(ligands)) {
			__ligand_idatm_types = ligands.get_idatm_types(__ligand_idatm_types);
			common::create_mols_from_seeds(__added, __seeds, ligands);
			ligands.clear();
		}
	}

	void FragmentLigands::__fragment_ligands( Molib::PDBreader& lpdb, const CmdLnOpts& cmdl ) {
		bool thread_is_not_done;
		Molib::Molecules ligands;
		{
			lock_guard<std::mutex> guard(__prevent_re_read_mtx);
			thread_is_not_done = lpdb.parse_molecule(ligands);
		}
		while(thread_is_not_done) {
			// Compute properties, such as idatm atom types, fragments, seeds,
			// rotatable bonds etc.
			ligands.compute_idatm_type()
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
			{
				lock_guard<std::mutex> guard(__add_to_typing_mtx);
				ligands.compute_overlapping_rigid_segments(cmdl.seeds_file());
				
				__ligand_idatm_types = ligands.get_idatm_types(__ligand_idatm_types);
				common::create_mols_from_seeds(__added, __seeds, ligands);
			}

			inout::output_file(ligands, cmdl.prep_file(), ios_base::app);
			ligands.clear();

			lock_guard<std::mutex> guard(__prevent_re_read_mtx);
			thread_is_not_done = lpdb.parse_molecule(ligands);
		}
	}

	void FragmentLigands::__continue_from_prev( const CmdLnOpts& cmdl) {

		cout << "Fragmenting files in " << cmdl.ligand_file() << endl;

		Molib::PDBreader lpdb(cmdl.ligand_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

		std::vector<std::thread> threads;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back( std::thread([&,this] {__fragment_ligands(lpdb, cmdl);} ) );
		}
		for(auto& thread : threads) {
			thread.join();
		}
		dbgmsg(__seeds);

		/* Erase atom properties since they make the graph matching incorrect
		 *
		 */
		__seeds.erase_properties();

	}
}
