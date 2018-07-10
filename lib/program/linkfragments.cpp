#include "candock/program/linkfragments.hpp"

#include <boost/filesystem.hpp>

#include "candock/linker/linker.hpp"
#include "candock/modeler/modeler.hpp"
#include "candock/helper/path.hpp"
#include "candock/helper/grep.hpp"

#include "candock/fileout/fileout.hpp"

namespace candock {
namespace Program {

        bool LinkFragments::__can_read_from_files () {
                std::regex regex;
                regex.assign ("REMARK   5 MOLECULE (\\w*)");
                std::ifstream file (cmdl.get_string_option ("prep"));

                std::vector<std::string> all_names = Grep::search_stream (file, regex);

                auto prot_name = boost::filesystem::basename(__receptor.name());
                boost::filesystem::path p (prot_name);
                p /= cmdl.get_string_option ("docked_dir");

                for (auto molec : all_names) {
                        boost::filesystem::path p2 = p / (molec + ".pdb");

                        if (Inout::file_size (p2.string()) <= 0) {
                                return false;
                        }
                }

                log_note << "Linking for all molecules in " << cmdl.get_string_option ("prep") << " for " << __receptor.name() << " is complete, skipping." << std::endl;

                return true;
        }

        void LinkFragments::__read_from_files () {

                std::regex regex;
                regex.assign ("REMARK   5 MOLECULE (\\w*)");
                std::ifstream file (cmdl.get_string_option ("prep"));

                std::vector<std::string> all_names = Grep::search_stream (file, regex);

                boost::filesystem::path p (__receptor.name());
                p /= cmdl.get_string_option ("docked_dir");

                for (auto molec : all_names) {
                        boost::filesystem::path p2 = p / (molec + ".pdb");

                        parser::FileParser conf (p2.string(), parser::skip_atom | parser::first_model, 1);
                        conf.parse_molecule (__all_top_poses);
                        __all_top_poses.last().set_name (molec);
                }

        }

        void LinkFragments::__link_ligand (molib::Molecule &ligand) {

                auto prot_name = boost::filesystem::basename(__receptor.name());
                boost::filesystem::path p (prot_name);
                p = p / cmdl.get_string_option ("docked_dir") / (ligand.name() + ".pdb");

                OMMIface::ForceField ffcopy (__ffield);

                // if docking of one ligand fails, docking of others shall continue...
                try {
                        ffcopy.insert_topology (ligand);

                        dbgmsg ("LINKING LIGAND : " << endl << ligand);

                        if (ligand.first().first().get_rigid().size() < 1) {
                                throw Error ("No seeds to link");
                        }

                        /**
                          * Read top seeds for this ligand
                          */
                        molib::NRset top_seeds = __seeds_database.get_seeds (ligand, cmdl.get_double_option ("top_percent"));

                        ligand.erase_properties(); // required for graph matching
                        top_seeds.erase_properties(); // required for graph matching

                        molib::Molecule crystal_ligand(ligand);

                        /**
                         * Jiggle the coordinates by one-thousand'th of an Angstrom to avoid minimization failures
                         * with initial bonded relaxation failed errors
                        */
                        std::mt19937 rng;
                        std::random_device::result_type seed;
                        
                        if (cmdl.get_int_option("jiggle_seed") != -1 && cmdl.get_int_option("jiggle_seed") != -2) {
                                seed = static_cast<std::random_device::result_type>(cmdl.get_int_option("jiggle_seed"));
                        } else {
                                // Fix issue with getting default values....
                                seed = std::random_device()() >> 2;
                                log_step << "Seed for " << ligand.name() << " is " << seed << std::endl;
                        }

                        if (cmdl.get_int_option("jiggle_seed") != -2) {
                            rng.seed(seed);
                            top_seeds.jiggle(rng);
                        }

                        /* Init minization options and constants, including ligand and receptor topology
                         *
                         */
                        OMMIface::Modeler modeler (
                                ffcopy,
                                cmdl.get_string_option ("fftype"),
                                cmdl.get_double_option ("mini_tol"),
                                cmdl.get_int_option ("max_iter"),
                                false
                        );

                        /**
                         * Connect seeds with rotatable linkers, account for symmetry, optimize
                         * seeds with appendices, minimize partial conformations between linking.
                         *
                         */

                        //FIXME Why does this modify (remove?) the ligand coordinates
                        linker::Linker linker (modeler, __receptor, ligand, top_seeds, __gridrec, __score,
                                               cmdl.get_bool_option ("cuda"), cmdl.get_bool_option ("iterative"),
                                               cmdl.get_int_option ("cutoff"), cmdl.get_int_option ("spin"),
                                               cmdl.get_double_option ("tol_seed_dist"), cmdl.get_double_option ("lower_tol_seed_dist"),
                                               cmdl.get_double_option ("upper_tol_seed_dist"),
                                               cmdl.get_int_option ("max_possible_conf"),
                                               cmdl.get_int_option ("link_iter"),
                                               cmdl.get_double_option ("clash_coeff"), cmdl.get_double_option ("docked_clus_rad"),
                                               cmdl.get_double_option ("max_allow_energy"), cmdl.get_int_option ("max_num_possibles"),
                                               cmdl.get_int_option ("max_clique_size"), cmdl.get_int_option ("max_iter_final"), cmdl.get_int_option ("max_iter_pre"),
                                               cmdl.get_string_option("platform"), cmdl.get_string_option("precision"), cmdl.get_string_option("accelerators"));

                        molib::Atom::Graph cryst_graph =  molib::Atom::create_graph(crystal_ligand.get_atoms());

                        linker::DockedConformation::Vec docks = linker.link();
                        linker::DockedConformation::sort (docks);

                        int model = 0;

                        for (auto &docked : docks) {

                                docked.get_ligand().change_residue_name ("CAN");

                                double rmsd = std::nan ("");

                                if (cmdl.get_bool_option ("rmsd_crystal")) {

                                        docked.get_ligand().erase_properties();
                                        molib::Atom::Vec atoms = docked.get_ligand().get_atoms();

                                        int reenum = 0;
                                        for ( auto &atom : atoms ) {
                                                atom->set_atom_number(++reenum);
                                                atom->erase_properties();
                                        }

                                        molib::Atom::Graph docked_graph = molib::Atom::create_graph(atoms);

                                        rmsd = molib::Atom::compute_rmsd(cryst_graph, docked_graph);
                                }

                                // output docked molecule conformations
                                std::stringstream ss;
                                ss << "REMARK   4 SEED USED FOR " << ligand.name() << " IS " << seed << "\n";
                                fileout::print_complex_pdb (ss, docked.get_ligand(), docked.get_receptor(),
                                                        docked.get_energy(), docked.get_potential_energy(),
                                                        ++model, docked.get_max_clq_identity(), rmsd);
                                Inout::output_file (ss.str(), p.string(), std::ios_base::app);
                        }

                        __all_top_poses.add (new molib::Molecule (docks[0].get_ligand()));
                        ffcopy.erase_topology (ligand);
                } catch (std::exception &e) {
                        log_error << "Error: skipping ligand " << ligand.name() << " with " << __receptor.name() << " due to : " << e.what() << std::endl;
                        std::stringstream ss;
                        ligand.change_residue_name ("CAN");
                        ss << "REMARK  20 non-binder " << ligand.name() << " with " << __receptor.name() << " because " << e.what() << std::endl << ligand;
                        Inout::file_open_put_stream (p.string(), ss, std::ios_base::app);
                }
        }

        void LinkFragments::__continue_from_prev () {

                log_step << "Starting to dock the fragments into originally given ligands" << std::endl;

                if (! Inout::file_size (cmdl.get_string_option ("prep"))) {
                        log_warning << cmdl.get_string_option ("prep") << " is either blank or missing, no (initial) ligand docking will take place.";
                        return;
                }

                parser::FileParser lpdb2 (cmdl.get_string_option ("prep"), parser::all_models, 1);

                std::mutex concurrent_numbering;
                std::mutex additon_to_top_dock;
                std::vector<std::thread> threads;

                for (int i = 0; i < cmdl.ncpu(); ++i) {
                        threads.push_back (std::thread ([ &,this] {
                                molib::Molecules ligands;

                                while (lpdb2.parse_molecule (ligands)) {

                                        molib::Molecule &ligand = ligands.first();
                                        auto prot_name = boost::filesystem::basename(__receptor.name());
                                        boost::filesystem::path p (prot_name);
                                        p = p / cmdl.get_string_option ("docked_dir") / (ligand.name() + ".pdb");

                                        if (Inout::file_size (p.string()) > 0) {
                                                std::lock_guard<std::mutex> guard(additon_to_top_dock);
                                                log_note << ligand.name() << " is alread docked to " << __receptor.name() << ", skipping." << std::endl;

                                                parser::FileParser conf (p.string(), parser::skip_atom | parser::first_model | parser::docked_poses_only, 1);
                                                conf.parse_molecule (__all_top_poses);
                                                __all_top_poses.last().set_name (ligand.name());
                                                ligands.clear();

                                                continue;
                                        }

                                        /**
                                         * Ligand's resn MUST BE UNIQUE for ffield
                                         */
                                        ligand.change_residue_name (concurrent_numbering, __ligand_cnt);

                                        __link_ligand (ligand);
                                        ligands.clear();
                                }
                        }));
                }

                for (auto &thread : threads) {
                        thread.join();
                }

                log_step << "Linking of fragments is complete" << std::endl;
        }

        void LinkFragments::link_ligands (const molib::Molecules &ligands) {
                size_t j = 0;

                std::vector<std::thread> threads;
                std::mutex counter_protect;
                std::mutex concurrent_numbering;
                std::mutex additon_to_top_dock;

                for (int i = 0; i < cmdl.ncpu(); ++i) {
                        threads.push_back (std::thread ([ &,this] {
                                while (true) {
                                        std::unique_lock<std::mutex> guard (counter_protect, std::defer_lock);
                                        guard.lock();

                                        if (j >= ligands.size())
                                                return;

                                        molib::Molecule &ligand = ligands[j++];

                                        guard.unlock();

                                        boost::filesystem::path p (__receptor.name());
                                        p = p / cmdl.get_string_option ("docked_dir") / (ligand.name() + ".pdb");

                                        if (Inout::file_size (p.string()) > 0) {
                                                std::lock_guard<std::mutex> guard(additon_to_top_dock);
                                                log_note << ligand.name() << " is alread docked to " << __receptor.name() << ", skipping." << std::endl;

                                                parser::FileParser conf (p.string(), parser::skip_atom | parser::first_model, 1);
                                                conf.parse_molecule (__all_top_poses);
                                                __all_top_poses.last().set_name (ligand.name());

                                                continue;
                                        }

                                        /**
                                         * Ligand's resn MUST BE UNIQUE for ffield
                                         */
                                        ligand.change_residue_name (concurrent_numbering, __ligand_cnt);
                                        __link_ligand (ligand);
                                }
                        }));
                }

                for (auto &thread : threads) {
                        thread.join();
                }
        }
}
}
