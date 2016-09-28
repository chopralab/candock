#include "target.hpp"

#include "pdbreader/pdbreader.hpp"

#include "findcentroids.hpp"

namespace Program {
	Target::Target(const std::string& input_name ) {

		/* Initialize parsers for receptor and read
		 * the receptor molecule(s)
		 * 
		 */
		for ( const auto &a : inout::Inout::files_matching_pattern (input_name, ".pdb")) {
			/* Initialize parsers for receptor (and ligands) and read
			 * the receptor molecule(s)
			 * 
			 */
			Molib::PDBreader rpdb(a, Molib::PDBreader::first_model);
			Molib::Molecules receptors = rpdb.parse_molecule();
			Molib::Molecule& current = __receptors.add(new Molib::Molecule ( std::move (receptors[0]) ));
			current.set_name( a.substr(0, a.length() - 4 ) );

			__preprecs.push_back(DockedReceptor {current, nullptr, nullptr, nullptr});
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
		for ( auto &a : __preprecs ) {
			std::unique_ptr<Molib::Atom::Grid> pgridrec(new Molib::Atom::Grid(a.protein.get_atoms()));
			a.gridrec = std::move(pgridrec);
		}
	}

	void Target::find_centroids(const CmdLnOpts& cmdl ) {
		for ( auto &a : __preprecs ) {
			std::unique_ptr<FindCentroids> pcen (new FindCentroids(a.protein));

			/* Run section of Candock designed to find binding site1s
			 * Currently, this runs ProBIS and does not require any
			 * previous step to be competed.
			 *
			 */

			pcen->run_step(cmdl);
			__preprecs.back().centroids = std::move(pcen);
		}
	}

	void Target::dock_fragments(const Molib::Score& score, const FragmentLigands& ligand_fragments, const CmdLnOpts& cmdl) {
		for ( auto &a : __preprecs ) {
			std::unique_ptr<DockFragments> pdockfragments(new DockFragments(*a.centroids, ligand_fragments, score, *a.gridrec, a.protein.name(), cmdl));
			pdockfragments->run_step(cmdl);
			a.prepseeds = std::move(pdockfragments);
		}
	}

	void Target::link_fragments(const Molib::Score& score, const Program::CmdLnOpts& cmdl) {
		__ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file())
			.parse_forcefield_file(cmdl.water_xml_file());

		for ( auto &a : __preprecs ) {
			a.protein.prepare_for_mm(__ffield, *a.gridrec);
			__ffield.insert_topology(a.protein);
			
			std::unique_ptr<LinkFragments> plinkfragments(new LinkFragments(a.protein, score, __ffield, *a.gridrec));
			plinkfragments->run_step(cmdl);
			a.dockedlig = std::move(plinkfragments);
		}
	}

	void Target::determine_overlapping_seeds(const int max_seeds, const int number_of_occurances, const bool is_anti_target) { //TODO: Remove cmdl dependency

		std::multiset<std::string> good_seed_list;

		for ( auto &a : __preprecs ) {
			auto result = a.prepseeds->get_best_seeds();

			if ( max_seeds != -1 && max_seeds < result.size()) 
				result.resize( max_seeds );

			for ( auto &b : result )
				good_seed_list.insert(b.second);
		}

		for ( auto c = good_seed_list.cbegin(); c != good_seed_list.cend(); ) {
			if ((static_cast<int>(good_seed_list.count(*c)) < number_of_occurances) && !is_anti_target) {
				c = good_seed_list.erase(c);
			} else {
				++c;
			}
		}
		
	}
}
