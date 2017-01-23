#include "target.hpp"

#include "pdbreader/pdbreader.hpp"

#include "helper/path.hpp"
#include "common.hpp"

namespace Program {
	Target::Target(const std::string& input_name ) {

		// If the user doesn't want to use this feature
		if (input_name == "")
			return;

		if (!boost::filesystem::exists(input_name)) {
			throw Error("Provided file or directory does not exist: " + input_name);
		}

		/* Initialize parsers for receptor and read
		 * the receptor molecule(s)
		 * 
		 */
		if ( inout::Inout::file_size(input_name) > 0 ) {
			// If the option given is a regular file, act like previous versions
			Molib::PDBreader rpdb(input_name, Molib::PDBreader::first_model);
			Molib::Molecules receptors = rpdb.parse_molecule();
			Molib::Molecule& current = __receptors.add(new Molib::Molecule ( std::move (receptors[0]) ));
			current.set_name(boost::filesystem::basename(input_name.substr(0, input_name.length() - 4))); // Emulate the original version of candock
                        boost::filesystem::create_directory(current.name());
			__preprecs.push_back(DockedReceptor {current, nullptr, nullptr, nullptr});
		} else for ( const auto &a : inout::Inout::files_matching_pattern (input_name, ".pdb")) {
			// Otherwise we treat it like the new version intends.
			Molib::PDBreader rpdb(a, Molib::PDBreader::first_model);
			Molib::Molecules receptors = rpdb.parse_molecule();
			Molib::Molecule& current = __receptors.add(new Molib::Molecule ( std::move (receptors[0]) ));
			current.set_name( a.substr(0, a.length() - 4 ) );
			boost::filesystem::create_directory(current.name());

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
			a.centroids = std::move(pcen);
		}
	}

	void Target::dock_fragments(const FragmentLigands& ligand_fragments, const CmdLnOpts& cmdl) {
		for ( auto &a : __preprecs ) {

			/* Read distributions file and initialize scores
			* 
			*/

			std::unique_ptr<Molib::Score> score (new Molib::Score(cmdl.ref_state(), cmdl.comp(),
																  cmdl.rad_or_raw(), cmdl.dist_cutoff(),
																  cmdl.step_non_bond()));

			score->define_composition(__receptors.get_idatm_types(),
									ligand_fragments.ligand_idatm_types())
									.process_distributions_file(cmdl.distributions_file())
									.compile_scoring_function()
									.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

			a.score = std::move(score);

			// Prepare the receptor for docking to
			std::unique_ptr<OMMIface::ForceField> ffield ( new OMMIface::ForceField);

			ffield->parse_gaff_dat_file(cmdl.gaff_dat_file())
				.add_kb_forcefield(*a.score, cmdl.step_non_bond())
				.parse_forcefield_file(cmdl.amber_xml_file())
				.parse_forcefield_file(cmdl.water_xml_file());

			a.ffield = std::move(ffield);

			a.protein.prepare_for_mm(*a.ffield, *a.gridrec);

			std::unique_ptr<DockFragments> pdockfragments(new DockFragments(*a.centroids, ligand_fragments, *a.score, *a.gridrec, a.protein.name(), cmdl));
			pdockfragments->run_step(cmdl);
			a.prepseeds = std::move(pdockfragments);
		}
	}

	void Target::link_fragments(const CmdLnOpts& cmdl) {

		for ( auto &a : __preprecs ) {
			a.ffield->insert_topology(a.protein);
			std::unique_ptr<LinkFragments> plinkfragments(new LinkFragments(a.protein, *a.score, *a.ffield, *a.gridrec));
			plinkfragments->run_step(cmdl);
			a.dockedlig = std::move(plinkfragments);
		}
	}

	void Target::design_ligands(const CmdLnOpts& cmdl, FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add ) {
#ifndef NDEBUG
		for (auto &s : seeds_to_add ) {
			cout << s << endl;
		}
#endif
		for ( auto &a : __preprecs ) {
			
			if ( inout::Inout::file_size("designed.pdb") ) {
				
				cout << "designed.pdb found -- skipping generation of new designs" << endl;
				
				Molib::PDBreader dpdb ("designed.pdb", Molib::PDBreader::all_models );
				Molib::Molecules designs;
				dpdb.parse_molecule(designs);

				ligand_fragments.add_seeds_from_molecules(designs, cmdl);
				a.prepseeds->run_step(cmdl);
				a.dockedlig->link_ligands(designs, cmdl);

				continue;
			}

			std::unique_ptr<design::Design> designer( new design::Design( a.dockedlig->top_poses().first() ));
			if (! seeds_to_add.empty() )
				designer->functionalize_hydrogens_with_fragments(common::read_top_seeds_files(seeds_to_add,
																 Path::join(a.protein.name(), cmdl.top_seeds_dir()),
																 cmdl.top_seeds_file(), cmdl.top_percent() ),
																 cmdl.tol_seed_dist(), cmdl.clash_coeff() );

			const vector<string>& h_single_atoms = cmdl.get_string_vector("add_single_atoms");
			const vector<string>& a_single_atoms = cmdl.get_string_vector("change_terminal_atom");

			if (!a_single_atoms.empty())
				designer->functionalize_extremes_with_single_atoms(a_single_atoms);
			if (!h_single_atoms.empty())
				designer->functionalize_hydrogens_with_single_atoms(h_single_atoms);
#ifndef NDEBUG
			inout::output_file(designer->get_internal_designs(), "internal_designs.pdb");
#endif
			inout::output_file(designer->prepare_designs(cmdl.seeds_file()), "designed.pdb");
			ligand_fragments.add_seeds_from_molecules(designer->designs(), cmdl);
			a.prepseeds->run_step(cmdl);
			a.dockedlig->link_ligands(designer->designs(), cmdl);
		}
	}

	std::multiset<std::string> Target::determine_overlapping_seeds(const int max_seeds, const int number_of_occurances) const {

		std::multiset<std::string> good_seed_list;

		for ( auto &a : __preprecs ) {
			auto result = a.prepseeds->get_best_seeds();

			if ( max_seeds != -1 && static_cast<size_t>(max_seeds) < result.size()) 
				result.resize( max_seeds );

			for ( auto &b : result )
				good_seed_list.insert(b.second);
		}

		for ( auto c = good_seed_list.cbegin(); c != good_seed_list.cend(); ) {
			if (static_cast<int>(good_seed_list.count(*c)) < number_of_occurances) {
				c = good_seed_list.erase(c);
			} else {
				++c;
			}
		}

		return good_seed_list;
	}

	std::set<std::string> Target::determine_non_overlapping_seeds( const Target& targets, const Target& antitargets, const CmdLnOpts& cmdl ) {
		set<string> solo_target_seeds;
		const vector<string>& forced_seeds = cmdl.get_string_vector("force_seed");

		if (forced_seeds.size() != 0 && forced_seeds[0] != "off") {
			std::copy( forced_seeds.begin(), forced_seeds.end(), std::inserter(solo_target_seeds, solo_target_seeds.end()));
		} else {
			cout << "Determining the best seeds to add" << endl;
			multiset<string>  target_seeds =     targets.determine_overlapping_seeds(cmdl.get_int_option("seeds_to_add"),   cmdl.get_int_option("seeds_till_good"));
			multiset<string> atarget_seeds = antitargets.determine_overlapping_seeds(cmdl.get_int_option("seeds_to_avoid"), cmdl.get_int_option("seeds_till_bad"));

			std::set_difference( target_seeds.begin(),  target_seeds.end(),
							    atarget_seeds.begin(), atarget_seeds.end(),
							    std::inserter(solo_target_seeds, solo_target_seeds.end())
			);
		}

		return solo_target_seeds;
	}
}
