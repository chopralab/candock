#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "program/cmdlnopts.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "program/common.hpp"
#include "helper/error.hpp"
#include "molib/grid.hpp"
#include "molib/nrset.hpp"
#include "molib/pdbreader.hpp"
#include "score/score.hpp"
#include "linker/linker.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"

#include "modeler/forcefield.hpp"
#include "score/score.hpp"
#include "modeler/modeler.hpp"

using namespace std;
using namespace Program;

////////////////// TEST OF RECEPTOR WITH COFACTORS ///////////////////////////

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv, CmdLnOpts::STARTING | CmdLnOpts::LIG_FRAMGENT);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		/* Initialize parser for receptor and read also cofactors
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

		/* Compute idatm atom types for receptor and cofactors
		 * 
		 */
		receptors.compute_idatm_type()
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

		inout::output_file(receptors, "prepared_receptor.pdb");
		
		cmdl.display_time("finished");

	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
