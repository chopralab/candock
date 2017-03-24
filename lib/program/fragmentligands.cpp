#include "fragmentligands.hpp"

#include <boost/filesystem.hpp>
#include <mutex>

#include "helper/inout.hpp"
#include "common.hpp"
#include "pdbreader/molecules.hpp"
#include "helper/options.hpp"

namespace Program {

	bool FragmentLigands::__can_read_from_files() {
		return (inout::Inout::file_size( cmdl.get_string_option("seeds_pdb") ) > 0 
		     || inout::Inout::file_size( cmdl.get_string_option("prep") ) > 0 )
                     && inout::Inout::file_size( cmdl.get_string_option("seeds") ) > 0;
	}

	void FragmentLigands::__read_from_files() {

		if ( inout::Inout::file_size( cmdl.get_string_option("seeds_pdb") ) <= 0 ) {

			cout << "Could not read seeds from " << cmdl.get_string_option("seeds_pdb") << endl;
			cout << "Reading fragmented files from " << cmdl.get_string_option("prep") << endl;

			Molib::PDBreader lpdb(cmdl.get_string_option("prep"), Molib::PDBreader::all_models, cmdl.get_int_option("max_num_ligands"));

			Molib::Molecules ligands;
			while(lpdb.parse_molecule(ligands)) {
				__ligand_idatm_types = ligands.get_idatm_types(__ligand_idatm_types);
				common::create_mols_from_seeds(__added, __seeds, ligands);
				ligands.clear();
			}

			__seeds.erase_properties();

			inout::output_file(__seeds, cmdl.get_string_option("seeds_pdb"), ios_base::out);
		} else {

			cout << "Reading seeds from: " << cmdl.get_string_option("seeds_pdb") << endl;

			Molib::PDBreader lpdb(cmdl.get_string_option("seeds_pdb"), Molib::PDBreader::all_models);
			lpdb.parse_molecule(__seeds);
			__ligand_idatm_types = __seeds.get_idatm_types(__ligand_idatm_types);
			for (const auto& mol : __seeds) {
				__added.insert( stoi(mol.name()) );
            }
		}
	}

	void FragmentLigands::__fragment_ligands( Molib::PDBreader& lpdb, const bool write_out_for_linking, const bool no_rotatable_bond) {
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

			if ( no_rotatable_bond ) {
				for (auto atom : ligands.get_atoms())
				for (auto bond : atom->get_bonds())
					bond->set_rotatable(""); // We still want other properties assigned with rotatable bonds to be assigned
			}

			lock_guard<std::mutex> guard(__add_to_typing_mtx);
			ligands.compute_overlapping_rigid_segments(cmdl.get_string_option("seeds"));

			__ligand_idatm_types = ligands.get_idatm_types(__ligand_idatm_types);
			common::create_mols_from_seeds(__added, __seeds, ligands);

			if (write_out_for_linking)
				inout::output_file(ligands, cmdl.get_string_option("prep"), ios_base::app);

			ligands.clear();
			thread_is_not_done = lpdb.parse_molecule(ligands);
		}
	}

	void FragmentLigands::__continue_from_prev() {

		if (inout::Inout::file_size(cmdl.get_string_option("ligand")) > 0) {
			cout << "Fragmenting files in " << cmdl.get_string_option("ligand") << endl;

			Molib::PDBreader lpdb(cmdl.get_string_option("ligand"), 
				Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
				cmdl.get_int_option("max_num_ligands"));

			std::vector<std::thread> threads1;

			for(int i = 0; i < cmdl.ncpu(); ++i) {
				threads1.push_back( std::thread([&,this] {__fragment_ligands(lpdb, true, false);} ) );
			}
			for(auto& thread : threads1) {
				thread.join();
			}
		}

		const string& fragment_bag = cmdl.get_string_option("fragment_bag");

		if (inout::Inout::file_size(fragment_bag) > 0) {

			cout << "Adding fragments from " << fragment_bag << endl;

			Molib::PDBreader lpdb_additional(fragment_bag, 
				Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
				cmdl.get_int_option("max_num_ligands"));

			std::vector<std::thread> threads2;

			for(int i = 0; i < cmdl.ncpu(); ++i) {
				threads2.push_back( std::thread([&,this] {__fragment_ligands(lpdb_additional, false, false);} ) );
			}
			for(auto& thread : threads2) {
				thread.join();
			}
		}

		const string& molecular_fragments = cmdl.get_string_option("fragment_mol");

		if (inout::Inout::file_size(molecular_fragments) > 0 ) {
			cout << "Adding molecular fragments from " << molecular_fragments << endl;

			Molib::PDBreader lpdb_additional(molecular_fragments, 
				Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
				cmdl.get_int_option("max_num_ligands"));

			std::vector<std::thread> threads3;

			for(int i = 0; i < cmdl.ncpu(); ++i) {
				threads3.push_back( std::thread([&,this] {__fragment_ligands(lpdb_additional, false, true);} ) );
			}
			for(auto& thread : threads3) {
				thread.join();
			}
		}

		dbgmsg(__seeds);

		/* Erase atom properties since they make the graph matching incorrect
		 *
		 */
		__seeds.erase_properties();

		inout::output_file(__seeds, cmdl.get_string_option("seeds_pdb"), ios_base::out);
	}

	void FragmentLigands::add_seeds_from_molecules(const Molib::Molecules& molecules) {
		__ligand_idatm_types = molecules.get_idatm_types(__ligand_idatm_types);
		common::create_mols_from_seeds(__added, __seeds, molecules);
		__seeds.erase_properties();
		inout::output_file(__seeds, cmdl.get_string_option("seeds_pdb"), ios_base::out);
	}
}
