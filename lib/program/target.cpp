#include "target.hpp"

#include <boost/filesystem.hpp>

#include "helper/path.hpp"
#include "fragmenter/unique.hpp"
#include "molib/molecule.hpp"
#include "parser/fileparser.hpp"

#include "modeler/systemtopology.hpp"
#include "modeler/modeler.hpp"

#include "findcentroids.hpp"
#include "dockfragments.hpp"
#include "linkfragments.hpp"

namespace Program {

        Target::DockedReceptor::~DockedReceptor() {
                if (score)
                        delete score;

                if (ffield)
                        delete ffield;

                if (gridrec)
                        delete gridrec;

                if (centroids)
                        delete centroids;

                if (prepseeds)
                        delete prepseeds;

                if (dockedlig)
                        delete dockedlig;
        }

        Target::Target (const std::string &input_name) {

                // If the user doesn't want to use this feature
                if (input_name == "")
                        return;

                if (!boost::filesystem::exists (input_name)) {
                        throw Error ("Provided file or directory does not exist: " + input_name);
                }

                /* Initialize parsers for receptor and read
                 * the receptor molecule(s)
                 *
                 */
                if (Inout::file_size (input_name) > 0) {
                        // If the option given is a regular file, act like previous versions
                        Parser::FileParser rpdb (input_name, Parser::first_model);
                        Molib::Molecules receptors = rpdb.parse_molecule();
                        Molib::Molecule &current = __receptors.add (new Molib::Molecule (std::move (receptors[0])));
                        current.set_name (boost::filesystem::basename (input_name.substr (0, input_name.length() - 4))); // Emulate the original version of candock
                        boost::filesystem::create_directory (current.name());
                        __preprecs.push_back (DockedReceptor (current, input_name));
                } else for (const auto &a : Inout::files_matching_pattern (input_name, ".pdb")) {
                                // Otherwise we treat it like the new version intends.
                                Parser::FileParser rpdb (a, Parser::first_model);
                                Molib::Molecules receptors = rpdb.parse_molecule();
                                Molib::Molecule &current = __receptors.add (new Molib::Molecule (std::move (receptors[0])));
                                current.set_name (a.substr (0, a.length() - 4));
                                boost::filesystem::create_directory (current.name());

                                __preprecs.push_back (DockedReceptor (current, a));
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
                for (auto &a : __preprecs) {
                        a.gridrec = new Molib::Atom::Grid (a.protein.get_atoms());
                }
        }

        /* Read distributions file and initialize scores
         * 
         */

        void Target::__initialize_score(const FragmentLigands &ligand_fragments) {
                for ( auto& a : __preprecs ) {
                        if ( a.score != nullptr) {
                                continue;
                        }

                        a.score = new Molib::Score(cmdl.get_string_option("ref"), cmdl.get_string_option("comp"),
                                                   cmdl.get_string_option("func"),cmdl.get_int_option("cutoff"),
                                                   cmdl.get_double_option("step"));

                        a.score->define_composition(a.protein.get_idatm_types(),
                                                   ligand_fragments.ligand_idatm_types())
                              .process_distributions_file(cmdl.get_string_option("dist"))
                              .compile_scoring_function()
                              .parse_objective_function(cmdl.get_string_option("obj_dir"), cmdl.get_double_option("scale"));
                }
        }

        void Target::__initialize_ffield() {
                for (auto& a : __preprecs) {
                        if ( a.ffield != nullptr) {
                                continue;
                        }

                        // Prepare the receptor for docking to
                        OMMIface::ForceField* ffield = new OMMIface::ForceField;

                        ffield->parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
                                .add_kb_forcefield(*a.score, cmdl.get_double_option("step"))
                                .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                                .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                        a.ffield = ffield;

                        a.protein.prepare_for_mm(*a.ffield, *a.gridrec);
                }
        }

        std::set<int> Target::get_idatm_types(const std::set<int>& previous) {
                return __receptors.get_idatm_types(previous);
        }

        void Target::find_centroids( ) {
                for ( auto &a : __preprecs ) {
                        
                        if ( a.centroids != nullptr ) {
                                continue;
                        }
                        
                        /* Run section of Candock designed to find binding site1s
                         * Currently, this runs ProBIS and does not require any
                         * previous step to be competed.
                         *
                         */

                        a.centroids = new FindCentroids(a.protein, a.filename);
                        a.centroids->run_step();
                }
        }

        void Target::rescore_docked(const FragmentLigands& ligand_fragments) {

                __initialize_score(ligand_fragments);
                __initialize_ffield();

                for ( auto &a : __preprecs ) {

                        Parser::FileParser lpdb(cmdl.get_string_option("prep"), 
                        Parser::all_models|Parser::hydrogens, 
                        cmdl.get_int_option("max_num_ligands"));

                        Molib::Molecules ligands = lpdb.parse_molecule();

                        for ( auto &ligand : ligands ) {
                                const double energy = a.score->non_bonded_energy(*a.gridrec, ligand);

                                Inout::output_file(Molib::Molecule::print_complex(ligand, a.protein, energy), 
                                Path::join( cmdl.get_string_option("docked_dir"), ligand.name() + ".pdb"),  ios_base::app);

                                cout << "Energy is " << energy << " for " << ligand.name() << " in " << a.protein.name() << endl;
                        }
                }
        }

        void Target::dock_fragments(const FragmentLigands& ligand_fragments) {

                find_centroids();

                __initialize_score(ligand_fragments);
                __initialize_ffield();

                for ( auto &a : __preprecs ) {

                        if (a.prepseeds != nullptr) {
                                continue;
                        }

                        a.prepseeds = new DockFragments(*a.centroids, ligand_fragments, *a.score, *a.gridrec, a.protein.name());
                        a.prepseeds->run_step();
                }
        }

        void Target::link_fragments(const FragmentLigands &ligand_fragments) {

                dock_fragments(ligand_fragments);

                OMMIface::SystemTopology::loadPlugins();

                for ( auto &a : __preprecs ) {

                        if ( a.dockedlig != nullptr ) {
                                continue;
                        }

                        a.ffield->insert_topology(a.protein);
                        a.dockedlig = new LinkFragments(a.protein, *a.score, *a.ffield, *a.prepseeds, *a.gridrec);
                        a.dockedlig->run_step();
                }
        }

        void Target::minimize_force(const FragmentLigands &ligand_fragments) {

                OMMIface::SystemTopology::loadPlugins();

                __initialize_score(ligand_fragments);
                __initialize_ffield();

                for ( auto &a : __preprecs ) {

                        Parser::FileParser lpdb(cmdl.get_string_option("prep"), Parser::all_models|Parser::hydrogens, -1);

                        Molib::Molecules ligands = lpdb.parse_molecule();

                        for (auto &ligand : ligands) {
                                try {

                                        /**
                                         * Minimize system
                                         */

                                        a.ffield->insert_topology (ligand);

                                        OMMIface::Modeler modeler (*a.ffield, cmdl.get_string_option("fftype"), cmdl.get_int_option("cutoff"),
                                                  cmdl.get_double_option("mini_tol"), cmdl.get_int_option("max_iter"), cmdl.get_int_option("update_freq"), 
                                                  cmdl.get_double_option("pos_tol"), false, 2.0);

                                        modeler.add_topology (a.protein.get_atoms());
                                        modeler.add_topology (ligand.get_atoms());

                                        modeler.init_openmm();

                                        modeler.add_crds (a.protein.get_atoms(), a.protein.get_crds());
                                        modeler.add_crds (ligand.get_atoms(), ligand.get_crds());

                                        modeler.init_openmm_positions();

                                        cout << "Initial energy for " << a.protein.name() << " and " << ligand.name() 
                                             << " = " << a.score->non_bonded_energy (*a.gridrec, ligand) << endl;

                                        modeler.minimize_state (ligand, a.protein, *a.score);

                                        // init with minimized coordinates
                                        Molib::Molecule minimized_receptor (a.protein, modeler.get_state (a.protein.get_atoms()));
                                        Molib::Molecule minimized_ligand (ligand, modeler.get_state (ligand.get_atoms()));

                                        minimized_receptor.undo_mm_specific();

                                        Molib::Atom::Grid gridrec (minimized_receptor.get_atoms());
                                        const double energy = a.score->non_bonded_energy (gridrec, minimized_ligand);

                                        cout << "Minimized energy = " << energy << endl;

                                        Inout::output_file (Molib::Molecule::print_complex (minimized_ligand, minimized_receptor, energy),
                                                            Path::join (cmdl.get_string_option("docked_dir"), minimized_ligand.name() + ".pdb"), ios_base::app);

                                } catch (exception &e) {
                                        log_error << "MINIMIZATION FAILED FOR LIGAND " << ligand.name()
                                                  << " because of " << e.what() << endl;
                                }

                                a.ffield->erase_topology (ligand);
                        }
                }
        }

        void Target::make_scaffolds(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add ) {
                Molib::Unique created_design("designed.txt");
                Molib::Molecules all_designs;

                log_step << "Starting iteration #0 (making a scaffold)" << endl;
                const string design_file = "designed_0.pdb";
                if ( Inout::file_size(design_file) ) {
                        log_note << design_file << " found -- skipping generation of new designs this iteration" << endl;
                        Parser::FileParser dpdb (design_file, Parser::all_models );
                        Molib::Molecules designs;
                        dpdb.parse_molecule(designs);

                        ligand_fragments.add_seeds_from_molecules(designs);
                        all_designs.add(designs);
                } else {
                    for (auto &a : __preprecs) {
                        Molib::NRset nr = a.prepseeds->get_top_seeds(seeds_to_add, cmdl.get_double_option("top_percent") );
                        for ( auto &molecules : nr ) {
                                design::Design designer (molecules.first(), created_design);
                                designer.change_original_name(molecules.name());
                                designer.functionalize_hydrogens_with_fragments(nr, cmdl.get_double_option("tol_seed_dist"), cmdl.get_double_option("clash_coeff"));
#ifndef NDEBUG
                Inout::output_file(designer.get_internal_designs(), "internal_designs.pdb", ios_base::app);
#endif
                                all_designs.add( designer.prepare_designs() );
                        }
                    }
                }

                if ( all_designs.size() == 0 ) {
                        log_step << "No new designs, exiting" << endl;
                        return;
                }

                all_designs .compute_hydrogen()
                            .compute_bond_order()
                            .compute_bond_gaff_type()
                            .refine_idatm_type()
                            .erase_hydrogen()  // needed because refine changes connectivities
                            .compute_hydrogen()   // needed because refine changes connectivities
                            .compute_ring_type()
                            .compute_gaff_type()
                            .compute_rotatable_bonds() // relies on hydrogens being assigned
                            .erase_hydrogen()
                            .compute_overlapping_rigid_segments(cmdl.get_string_option("seeds"));

                Inout::output_file(all_designs, design_file);

                ligand_fragments.add_seeds_from_molecules(all_designs);
                for ( auto &b : __preprecs ) {
                        b.prepseeds->run_step();
                        b.dockedlig->clear_top_poses();
                        b.dockedlig->link_ligands(all_designs);
                }
        }

        void Target::design_ligands(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add ) {
#ifndef NDEBUG
                for (auto &s : seeds_to_add ) {
                        cout << s << endl;
                }
#endif

                int n = 1;
                Molib::Unique created_design("designed.txt");
                while (true) { // Probably a bad idea

                    log_step << "Starting design iteration #" << n << endl;

                    string design_file = "designed_" + std::to_string(n++) + ".pdb";
                    Molib::Molecules all_designs;
                    if ( Inout::file_size(design_file) ) {
                            log_note << design_file << " found -- skipping generation of new designs this iteration" << endl;
                            Parser::FileParser dpdb (design_file, Parser::all_models );
                            Molib::Molecules designs;
                            dpdb.parse_molecule(designs);

                            ligand_fragments.add_seeds_from_molecules(designs);
                            all_designs.add(designs);
                    } else {
                        for ( auto &a : __preprecs ) {
                            for ( auto &molecule : a.dockedlig->top_poses() ) {
                                design::Design designer ( molecule, created_design);
                                if (! seeds_to_add.empty() )
                                        designer.functionalize_hydrogens_with_fragments(a.prepseeds->get_top_seeds(seeds_to_add,cmdl.get_double_option("top_percent")),
                                                                                        cmdl.get_double_option("tol_seed_dist"), cmdl.get_double_option("clash_coeff") );

                                const vector<string>& h_single_atoms = cmdl.get_string_vector("add_single_atoms");
                                const vector<string>& a_single_atoms = cmdl.get_string_vector("change_terminal_atom");

                                if (!a_single_atoms.empty())
                                        designer.functionalize_extremes_with_single_atoms(a_single_atoms);
                                if (!h_single_atoms.empty())
                                        designer.functionalize_hydrogens_with_single_atoms(h_single_atoms);
#ifndef NDEBUG
                                Inout::output_file(designer.get_internal_designs(), "internal_designs.pdb", ios_base::app);
#endif
                                all_designs.add( designer.prepare_designs() );
                            }
                        }
                    }

                    if ( all_designs.size() == 0 ) {
                            log_step << "No new designs, exiting" << endl;
                            return;
                    }

                    all_designs.compute_hydrogen()
                               .compute_bond_order()
                               .compute_bond_gaff_type()
                               .refine_idatm_type()
                               .erase_hydrogen()  // needed because refine changes connectivities
                               .compute_hydrogen()   // needed because refine changes connectivities
                               .compute_ring_type()
                               .compute_gaff_type()
                               .compute_rotatable_bonds() // relies on hydrogens being assigned
                               .erase_hydrogen()
                               .compute_overlapping_rigid_segments(cmdl.get_string_option("seeds"));

                    Inout::output_file(all_designs, design_file);

                    ligand_fragments.add_seeds_from_molecules(all_designs);
                    for ( auto &b : __preprecs ) {
                        b.prepseeds->run_step();
                        b.dockedlig->clear_top_poses();
                        b.dockedlig->link_ligands(all_designs);
                    }
                }
        }

        void Target::make_objective() {
                Molib::Score score(cmdl.get_string_option("ref"), "complete",
                                   cmdl.get_string_option("func"), cmdl.get_int_option("cutoff"),
                                   cmdl.get_double_option("step"));

                score.define_composition(set<int>(), set<int>())
                     .process_distributions_file(cmdl.get_string_option("dist"))
                     .compile_objective_function();

                score.output_objective_function(cmdl.get_string_option("obj_dir"));

                Inout::output_file(score, cmdl.get_string_option("potential_file"));

        }

        std::multiset<std::string> Target::determine_overlapping_seeds (const int max_seeds, const int number_of_occurances) const {

                std::multiset<std::string> good_seed_list;

                for (auto &a : __preprecs) {
                        auto result = a.prepseeds->get_best_seeds();

                        if (max_seeds != -1 && static_cast<size_t> (max_seeds) < result.size())
                                result.resize (max_seeds);

                        for (auto &b : result)
                                good_seed_list.insert (b.second);
                }

                for (auto c = good_seed_list.cbegin(); c != good_seed_list.cend();) {
                        if (static_cast<int> (good_seed_list.count (*c)) < number_of_occurances) {
                                c = good_seed_list.erase (c);
                        } else {
                                ++c;
                        }
                }

                return good_seed_list;
        }

        std::set<std::string> Target::determine_non_overlapping_seeds (const Target &targets, const Target &antitargets) {
                set<string> solo_target_seeds;
                const vector<string> &forced_seeds = cmdl.get_string_vector ("force_seed");

                if (forced_seeds.size() != 0 && forced_seeds[0] != "off") {
                        std::copy (forced_seeds.begin(), forced_seeds.end(), std::inserter (solo_target_seeds, solo_target_seeds.end()));
                } else {
                        log_step << "Determining the best seeds to add" << endl;
                        multiset<string>  target_seeds =     targets.determine_overlapping_seeds (cmdl.get_int_option ("seeds_to_add"),   cmdl.get_int_option ("seeds_till_good"));
                        multiset<string> atarget_seeds = antitargets.determine_overlapping_seeds (cmdl.get_int_option ("seeds_to_avoid"), cmdl.get_int_option ("seeds_till_bad"));

                        std::set_difference (target_seeds.begin(),  target_seeds.end(),
                                             atarget_seeds.begin(), atarget_seeds.end(),
                                             std::inserter (solo_target_seeds, solo_target_seeds.end())
                                            );
                }

                return solo_target_seeds;
        }
}
