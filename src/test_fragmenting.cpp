#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "opts_candock.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "helper/error.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecules.hpp"
#include "pdbreader/pdbreader.hpp"

using namespace std;
using namespace Program;

////////////////// TEST FRAGMENTING OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		Molib::PDBreader lpdb(cmdl.ligand_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

		vector<thread> threads;
		mutex mtx;
		mutex mtx_read_lock;

		Molib::Molecules seeds;
		set<int> added;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb, &seeds, &added, &mtx, &mtx_read_lock, &cmdl] () {
				bool thread_is_not_done;
				Molib::Molecules ligands;
				{
					lock_guard<std::mutex> guard(mtx_read_lock);
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
						lock_guard<std::mutex> guard(mtx);
						ligands.compute_overlapping_rigid_segments(cmdl.seeds_file());
					}

					inout::output_file(ligands, cmdl.prep_file(), ios_base::app);
					ligands.clear();

					lock_guard<std::mutex> guard(mtx_read_lock);
					thread_is_not_done = lpdb.parse_molecule(ligands);
				}
			}));
		}
		for(auto& thread : threads) {
			thread.join();
		}
		
		dbgmsg("SEEDS FOUND IN THE MOLECULE : " << endl << seeds);
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
