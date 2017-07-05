#include "linkfragments.hpp"

#include <boost/filesystem.hpp>

#include "linker/linker.hpp"
#include "modeler/modeler.hpp"
#include "helper/path.hpp"
#include "helper/grep.hpp"

#include "fileout/fileout.hpp"

namespace Program {

        bool LinkFragments::__can_read_from_files () {
                boost::regex regex;
                regex.assign ("REMARK   5 MOLECULE (\\w*)");
                std::ifstream file (cmdl.get_string_option ("prep"));

                std::vector<std::string> all_names = Grep::search_stream (file, regex);

                boost::filesystem::path p (__receptor.name());
                p /= cmdl.get_string_option ("docked_dir");

                bool is_done = true;

                for (auto molec : all_names) {
                        boost::filesystem::path p2 = p / (molec + ".pdb");

                        if (Inout::file_size (p2.string()) <= 0) {
                                is_done = false;
                        } else {
                                Parser::FileParser conf (p2.string(), Parser::skip_atom | Parser::first_model, 1);
                                conf.parse_molecule (__all_top_poses);
                                __all_top_poses.last().set_name (molec);
                        }
                }

                return is_done;
        }

        void LinkFragments::__read_from_files () {
                log_note << "Linking for all molecules in " << cmdl.get_string_option ("prep") << " for " << __receptor.name() << " is complete, skipping." << endl;
        }

        void LinkFragments::__link_ligand (Molib::Molecule &ligand, const OMMIface::ForceField &ffield) {
                boost::filesystem::path p (__receptor.name());
                p = p / cmdl.get_string_option ("docked_dir") / (ligand.name() + ".pdb");

                if (Inout::file_size (p.string()) > 0) {
                        log_note << ligand.name() << " is alread docked to " << __receptor.name() << ", skipping." << endl;
                        return;
                }

                // if docking of one ligand fails, docking of others shall continue...
                try {
                        dbgmsg ("LINKING LIGAND : " << endl << ligand);

                        if (ligand.first().first().get_rigid().size() < 1) {
                                throw Error ("No seeds to link");
                        }

                        /**
                          * Read top seeds for this ligand
                          */
                        Molib::NRset top_seeds = __seeds_database.get_top_seeds (ligand, cmdl.get_double_option ("top_percent"));

                        ligand.erase_properties(); // required for graph matching
                        top_seeds.erase_properties(); // required for graph matching

                        Molib::Molecule crystal_ligand(ligand);

                        /**
                         * Jiggle the coordinates by one-thousand'th of an Angstrom to avoid minimization failures
                         * with initial bonded relaxation failed errors
                        */
                        top_seeds.jiggle();

                        /* Init minization options and constants, including ligand and receptor topology
                         *
                         */
                        OMMIface::Modeler modeler (ffield, cmdl.get_string_option ("fftype"), cmdl.get_int_option ("cutoff"),
                                                   cmdl.get_double_option ("mini_tol"), cmdl.get_int_option ("max_iter"), cmdl.get_int_option ("update_freq"),
                                                   cmdl.get_double_option ("pos_tol"), false, cmdl.get_double_option("dynamic_step_size"),
                                                   cmdl.get_double_option("temperature"), cmdl.get_double_option("friction"));

                        /**
                         * Connect seeds with rotatable linkers, account for symmetry, optimize
                         * seeds with appendices, minimize partial conformations between linking.
                         *
                         */

                        //FIXME Why does this modify (remove?) the ligand coordinates
                        Linker::Linker linker (modeler, __receptor, ligand, top_seeds, __gridrec, __score,
                                               cmdl.get_bool_option ("cuda"), cmdl.get_bool_option ("iterative"),
                                               cmdl.get_int_option ("cutoff"), cmdl.get_int_option ("spin"),
                                               cmdl.get_double_option ("tol_seed_dist"), cmdl.get_double_option ("lower_tol_seed_dist"),
                                               cmdl.get_double_option ("upper_tol_seed_dist"),
                                               cmdl.get_int_option ("max_possible_conf"),
                                               cmdl.get_int_option ("link_iter"),
                                               cmdl.get_double_option ("clash_coeff"), cmdl.get_double_option ("docked_clus_rad"),
                                               cmdl.get_double_option ("max_allow_energy"), cmdl.get_int_option ("max_num_possibles"),
                                               cmdl.get_int_option ("max_clique_size"), cmdl.get_int_option ("max_iter_final"));

                        Molib::Atom::Graph cryst_graph =  Molib::Atom::create_graph(crystal_ligand.get_atoms());

                        Linker::DockedConformation::Vec docks = linker.link();
                        Linker::DockedConformation::sort (docks);

                        int model = 0;

                        for (auto &docked : docks) {

                                docked.get_ligand().change_residue_name ("CAN");

                                double rmsd = std::nan ("");

                                if (cmdl.get_bool_option ("rmsd_crystal")) {

                                        docked.get_ligand().erase_properties();
                                        Molib::Atom::Vec atoms = docked.get_ligand().get_atoms();

                                        int reenum = 0;
                                        for ( auto &atom : atoms ) {
                                                atom->set_atom_number(++reenum);
                                                atom->erase_properties();
                                        }

                                        Molib::Atom::Graph docked_graph = Molib::Atom::create_graph(atoms);

                                        rmsd = Molib::Atom::compute_rmsd(cryst_graph, docked_graph);
                                }

                                // output docked molecule conformations
                                std::stringstream ss;
                                Fileout::print_complex_pdb (ss, docked.get_ligand(), docked.get_receptor(),
                                                        docked.get_energy(), docked.get_potential_energy(),
                                                        ++model, docked.get_max_clq_identity(), rmsd);
                                Inout::output_file (ss.str(), p.string(), ios_base::app);
                        }

                        __all_top_poses.add (new Molib::Molecule (docks[0].get_ligand()));
                } catch (exception &e) {
                        log_error << "Error: skipping ligand " << ligand.name() << " with " << __receptor.name() << " due to : " << e.what() << endl;
                        stringstream ss;
                        ligand.change_residue_name ("CAN");
                        ss << "REMARK  20 non-binder " << ligand.name() << " with " << __receptor.name() << " because " << e.what() << endl << ligand;
                        Inout::file_open_put_stream (p.string(), ss, ios_base::app);
                }
        }

        void LinkFragments::__continue_from_prev () {

                log_step << "Starting to dock the fragments into originally given ligands" << endl;

                if (! Inout::file_size (cmdl.get_string_option ("prep"))) {
                        log_warning << cmdl.get_string_option ("prep") << " is either blank or missing, no (initial) ligand docking will take place.";
                        return;
                }

                Parser::FileParser lpdb2 (cmdl.get_string_option ("prep"), Parser::all_models, 1);

                std::vector<std::thread> threads;

                for (int i = 0; i < cmdl.ncpu(); ++i) {
                        threads.push_back (std::thread ([ &,this] {
                                OMMIface::ForceField ffcopy (__ffield);
                                Molib::Molecules ligands;

                                while (lpdb2.parse_molecule (ligands)) {

                                        Molib::Molecule &ligand = ligands.first();
                                        /**
                                         * Ligand's resn MUST BE UNIQUE for ffield
                                         */
                                        ligand.change_residue_name (__concurrent_numbering, __ligand_cnt);

                                        try {
                                                ffcopy.insert_topology (ligand);
                                                __link_ligand (ligand, ffcopy);
                                                ffcopy.erase_topology (ligand); // he he
                                                ligands.clear();
                                        } catch (exception &e) {
                                                log_error << "Error: problem with ligand " << ligand.name() << " due to : " << e.what() << endl;
                                        }
                                }
                        }));
                }

                for (auto &thread : threads) {
                        thread.join();
                }

                log_step << "Linking of fragments is complete" << endl;
        }

        void LinkFragments::link_ligands (const Molib::Molecules &ligands) {
                size_t j = 0;

                std::vector<std::thread> threads;
                std::mutex counter_protect;

                for (int i = 0; i < cmdl.ncpu(); ++i) {
                        threads.push_back (std::thread ([ &,this] {
                                OMMIface::ForceField ffcopy (__ffield);

                                while (true) {
                                        unique_lock<std::mutex> guard (counter_protect, std::defer_lock);
                                        guard.lock();

                                        if (j >= ligands.size())
                                                return;

                                        Molib::Molecule &ligand = ligands[j++];

                                        guard.unlock();

                                        /**
                                         * Ligand's resn MUST BE UNIQUE for ffield
                                         */
                                        ligand.change_residue_name (__concurrent_numbering, __ligand_cnt);
                                        ffcopy.insert_topology (ligand);
                                        __link_ligand (ligand, ffcopy);
                                        ffcopy.erase_topology (ligand); // he he
                                }
                        }));
                }

                for (auto &thread : threads) {
                        thread.join();
                }
        }
}
