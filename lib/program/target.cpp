/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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

#include "fileout/fileout.hpp"

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
                              .compile_scoring_function();
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
                                .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                                .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                        if ( ! cmdl.get_string_option("gaff_heme").empty() ) {
                                dbgmsg( "Adding " << cmdl.get_string_option("gaff_heme") << endl);
                                ffield->parse_gaff_dat_file(cmdl.get_string_option("gaff_heme"));
                        }

                        a.ffield = ffield;

                        a.protein.prepare_for_mm(*a.ffield, *a.gridrec);
                }
        }

        void Target::__initialize_kbforce() {
                OMMIface::SystemTopology::loadPlugins();

                for (auto& a : __preprecs) {
                        a.score->parse_objective_function(cmdl.get_string_option("obj_dir"), cmdl.get_double_option("scale"), 1501);
                        a.ffield->add_kb_forcefield(*a.score, cmdl.get_double_option("step"), 15);
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

                __initialize_kbforce();

                for ( auto &a : __preprecs ) {

                        if ( a.dockedlig != nullptr ) {
                                continue;
                        }

                        a.ffield->insert_topology(a.protein);
                        a.dockedlig = new LinkFragments(a.protein, *a.score, *a.ffield, *a.prepseeds, *a.gridrec);
                        a.dockedlig->run_step();
                }
        }

        void Target::make_scaffolds(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add ) {
            
                std::stringstream used_seeds;
                for (auto &s : seeds_to_add ) {
                        used_seeds << s << endl;
                }
                
                Inout::output_file (used_seeds.str(),"new_scaffold_seeds.lst");

            
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

                created_design.write_out();

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
                
                std::stringstream ss;
                for (const auto &m : all_designs) {
                        Fileout::print_mol2(ss,m);
                }
                Inout::output_file(ss.str(), "designed_0.mol2");

                ligand_fragments.add_seeds_from_molecules(all_designs);
                for ( auto &b : __preprecs ) {
                        b.prepseeds->run_step();
                        b.dockedlig->clear_top_poses();
                        b.dockedlig->link_ligands(all_designs);
                }
        }

        void Target::design_ligands(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add ) {

                std::stringstream used_seeds;
                for (auto &s : seeds_to_add ) {
                        used_seeds << s << endl;
                }
                
                Inout::output_file (used_seeds.str(),"lead_optimization_seeds.lst");

                int n = 0;
                Molib::Unique created_design("designed.txt");
                while (true) { // Probably a bad idea

                    log_step << "Starting design iteration #" << ++n << endl;

                    string design_file = "designed_" + std::to_string(n) + ".pdb";
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

                    created_design.write_out();

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

                    std::stringstream ss;
                    for (const auto &m : all_designs) {
                            Fileout::print_mol2(ss,m);
                    }
                    Inout::output_file(ss.str(), "designed_" + std::to_string(n) + ".mol2");

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
