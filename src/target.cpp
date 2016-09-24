#include "target.hpp"

#include "pdbreader/pdbreader.hpp"

#include "findcentroids.hpp"

namespace Program {
	Target::Target(const CmdLnOpts& cmdl) {

		/* Initialize parsers for receptor and read
		 * the receptor molecule(s)
		 * 
		 */
		if (!cmdl.get_string_option("target_dir").empty()) {
			for ( const auto &a : inout::Inout::files_matching_pattern (cmdl.get_string_option("target_dir"), ".pdb")) {
				/* Initialize parsers for receptor (and ligands) and read
				 * the receptor molecule(s)
				 * 
				 */
				Molib::PDBreader rpdb(a, Molib::PDBreader::first_model);
				Molib::Molecules receptors = rpdb.parse_molecule();
				Molib::Molecule& current = __receptors.add(new Molib::Molecule (receptors[0]));
				current.set_name(a.substr(0, a.length() - 4 ) );
				cout << a.substr(0, a.length() - 4 ) << endl;

				/* Run section of Candock designed to find binding site1s
				 * Currently, this runs ProBIS and does not require any
				 * previous step to be competed.
				 *
				 */
				__centroids.push_back( FindCentroids(current) );
				__centroids.back().run_step(cmdl);
			}
		}

		/* Compute atom types for receptor and cofactor(s): gaff types for protein, 
		 * Mg ions, and water are read from the forcefield xml file later on while 
		 * gaff types for cofactors (ADP, POB, etc.) are calculated de-novo here
		 * 
		 */
		__receptors.compute_idatm_type()
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

		/* Create receptor grid
		 * 
		 */
		__gridrec = Molib::Atom::Grid(__receptors[0].get_atoms());

	}
}
