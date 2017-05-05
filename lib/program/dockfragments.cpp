#include "dockfragments.hpp"

#include <boost/filesystem/path.hpp>

#include "helper/path.hpp"
#include "helper/inout.hpp"
#include "options.hpp"
#include "helper/grep.hpp"

namespace Program {

        DockFragments::DockFragments( const FindCentroids& found_centroids,
                                      const FragmentLigands& fragmented_ligands,
                                      const Molib::Score& score,
                                      const Molib::Atom::Grid& gridrec,
                                      const std::string& name ) :
                                      __found_centroids(found_centroids),
                                      __fragmented_ligands(fragmented_ligands),
                                      __score(score),
                                      __gridrec(gridrec),
                                      __name(name) {
                if ( ! cmdl.get_string_option("top_seeds_dir").empty() ) {
                        __top_seeds_location = cmdl.get_string_option("top_seeds_dir");
                        return;
                }

                boost::filesystem::path p (__name);
                p = p / "top_seeds";
                __top_seeds_location = p.string();
        }

        bool DockFragments::__can_read_from_files () {
                const Molib::Molecules& all_seeds = __fragmented_ligands.seeds();

                // No early return so that we have the ability to redock missing seeds
                const string &top_seeds_file= cmdl.get_string_option("top_seeds_file");
                for (const auto& seed : all_seeds) {
                        boost::filesystem::path p (__top_seeds_location);
                        p = p / seed.name() / top_seeds_file;
                        if (Inout::file_size(p.string()) <= 0)
                                return false;
                }

                return true;
        }

        void DockFragments::__read_from_files () {
                log_step << "All seeds are present in " << cmdl.get_string_option("top_seeds_dir") << " for " << __name << ". Docking of fragments skipped." << endl;
        }

        void DockFragments::__dock_fragment ( int start, const Docker::Gpoints& gpoints, const Docker::Gpoints& gpoints0 ) {
                // iterate over docked seeds and dock unique seeds
                const string &top_seeds_file= cmdl.get_string_option("top_seeds_file");

                for (size_t j = start; j < __fragmented_ligands.seeds().size(); j+= cmdl.ncpu()) {
                        try {
                                boost::filesystem::path p (__top_seeds_location);
                                p = p / __fragmented_ligands.seeds()[j].name() / top_seeds_file;

                                if ( Inout::file_size(p.string()) > 0 ) {
                                        log_note << "Skipping docking of seed: " << __fragmented_ligands.seeds()[j].name() << " because it is already docked!" << endl;
                                        continue;
                                } else {
                                        log_note << "Docking seed: " << __fragmented_ligands.seeds()[j].name() << endl;
                                }
                                dbgmsg(__fragmented_ligands.seeds()[j]);
                                /* Compute all conformations of this seed with the center
                                 * atom fixed on coordinate origin using maximum clique algorithm
                                 * 
                                 */
                                Docker::Conformations conf(__fragmented_ligands.seeds()[j], gpoints0, 
                                                            cmdl.get_double_option("conf_spin"), // degrees
                                                            cmdl.get_int_option("num_univec") // number of unit vectors
                                                          );

#ifndef NDEBUG
                                Inout::output_file(conf, "conf_" + __fragmented_ligands.seeds()[j].name() + ".pdb"); 
#endif
                                /* Dock this seed's conformations to the entire grid by moving them 
                                 * over all gridpoints and probe where they clash with the receptor: 
                                 * cluster docked conformations based on rmsd and energy and output 
                                 * only best-scored cluster representatives
                                 *
                                 */
#ifndef NDEBUG
                                Docker::Dock dock(gpoints, conf, __fragmented_ligands.seeds()[j], __score, __gridrec,  cmdl.get_double_option("clus_rad"));
#else
                                Docker::Dock dock(gpoints, conf, __fragmented_ligands.seeds()[j], cmdl.get_double_option("clus_rad"));
#endif

                                dock.run();

                                Inout::output_file(dock.get_docked(), p.string()); // output docked & clustered seeds
                        }
                        catch (exception& e) {
                                cerr << "skipping seed due to : " << e.what() << endl;
                        }
                }
        }

        void DockFragments::__continue_from_prev () {

                log_step << "Docking fragments into: " << __top_seeds_location << 
                        ". Files will be named: " << cmdl.get_string_option("top_seeds_file") << endl;

                /* Create gridpoints for ALL centroids representing one or more binding sites
                 * 
                 */
                Docker::Gpoints gpoints(__score, __fragmented_ligands.ligand_idatm_types(), __found_centroids.centroids(),
                                        __gridrec, cmdl.get_double_option("grid"), cmdl.get_int_option("cutoff"),
                                        cmdl.get_double_option("excluded"), 
                                        cmdl.get_double_option("interatomic"));
                Inout::output_file(gpoints, Path::join(__name, cmdl.get_string_option("gridpdb_hcp")));

                /* Create a zero centered centroid with 10 A radius (max fragment 
                 * radius) for getting all conformations of each seed
                 * 
                 */
                Docker::Gpoints gpoints0(cmdl.get_double_option("grid"), cmdl.get_double_option("max_frag_radius"));

                /* Create template grids using ProBiS-ligands algorithm
                 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
                 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
                 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
                 */
                std::vector<std::thread> threads;

                for(int i = 0; i < cmdl.ncpu(); ++i) {
                        threads.push_back( std::thread([&,this,i] {__dock_fragment(i, gpoints, gpoints0);} ) );
                }
                for(auto& thread : threads) {
                        thread.join();
                }
                
                log_step << "Done with fragment docking" << std::endl;

        }

        std::vector<std::pair<double, std::string>> DockFragments::get_best_seeds() const {
                const Molib::Molecules& all_seeds = __fragmented_ligands.seeds();

                std::vector<std::pair<double, std::string>> seed_score_map;

                for ( const auto &seed : all_seeds ) {
                        //FIXME: Do not read from disk here
                        dbgmsg("Reading: " << seed.name() << endl);
                        try {
                                boost::filesystem::path p (__top_seeds_location);
                                p = p / seed.name() / cmdl.get_string_option("top_seeds_file");
                                Parser::FileParser spdb (p.string(), Parser::first_model, 1);

                                Molib::Molecules seed_molec = spdb.parse_molecule();
                                seed_score_map.push_back( {std::stod( seed_molec.first().name()), seed.name()} );
                        } catch ( Error e) {
                                        cerr << "Skipping seed " << seed.name() << " in " << __name << " because " << e.what() << endl;
                        }
                }

                std::sort(seed_score_map.begin(), seed_score_map.end() );
                return seed_score_map;
        }

        Molib::NRset DockFragments::get_top_seeds(const std::set<std::string> &seeds, const double top_percent) const {
                Molib::NRset top_seeds;
                
                boost::regex regex;
                regex.assign("REMARK   5 MOLECULE ", boost::regex_constants::basic);

                for ( auto &fragment : seeds ) {

                        boost::filesystem::path file_to_read(__top_seeds_location);
                        file_to_read /= fragment;
                        file_to_read /= cmdl.get_string_option("top_seeds_file");
                        std::ifstream file( file_to_read.c_str() );

                        const size_t number_of_seeds = Grep::count_matches(file, regex);
                        const int sz = static_cast<int> (number_of_seeds * top_percent);
                        dbgmsg("taking " << sz << " top seeds for seed " << fragment);

                        // Add one in case the user is silly enough to select a top_percent of 0.000
                        Parser::FileParser pdb( file_to_read.string(), Parser::all_models, sz + 1 );

                        dbgmsg("reading top_seeds_file for seed id = " << fragment);
                        Molib::Molecules all = pdb.parse_molecule();

                        dbgmsg("number of top seeds left = " << all.size());

                        Molib::Molecules &last = top_seeds.add(new Molib::Molecules(all));

                        if (last.empty()) {
                                throw Error("die : there are no docked conformations for seed " + fragment);
                        }
                }

                return top_seeds;
        }

        Molib::NRset DockFragments::get_top_seeds(const Molib::Molecule &ligand, const double top_percent) const {
                std::set<string> seeds_to_read;
                const Molib::Model &model = ligand.first().first();
                for (auto &fragment : model.get_rigid()) { // iterate over seeds
                        if (fragment.is_seed()) {
                                seeds_to_read.insert(std::to_string(fragment.get_seed_id()));
                        }
                }
                return get_top_seeds(seeds_to_read, top_percent);
        }
}
