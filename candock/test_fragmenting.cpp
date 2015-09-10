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

		vector<thread> threads;
		mutex mtx;

		Molib::Molecules seeds;
		set<int> added;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb, &seeds, &added, &mtx] () {
				Molib::Molecules ligands;
				while(lpdb.parse_molecule(ligands)) {
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
