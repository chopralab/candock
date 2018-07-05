#include "candock/program/cmdlnopts.hpp"
#include "candock/program/fragmentligands.hpp"
#include "candock/program/target.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"

#include "candock/score/score.hpp"
#include "candock/modeler/forcefield.hpp"

#include <memory>

using namespace std;
using namespace candock;
using namespace candock::Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
        try {

                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                help::Options::set_options(new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::FORCE_FIELD| 
                                Program::CmdLnOpts::SCORING));

                cerr << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                cerr << help::Options::get_options()->configuration_file() << endl;

                molib::Molecules receptor_mols;
                molib::Molecules ligand_mols;

                parser::FileParser drpdb (cmdl.get_string_option("receptor"),
                                          parser::pdb_read_options::protein_poses_only |
                                          parser::pdb_read_options::all_models);

                drpdb.parse_molecule(receptor_mols);
                receptor_mols.compute_idatm_type();
                receptor_mols.erase_hydrogen();

                // Check to see if the user set the ligand option, use the receptor as a last resort.
                const std::string &ligand_file = Inout::file_size(cmdl.get_string_option("ligand")) == 0 ?
                                                    cmdl.get_string_option("receptor") :
                                                    cmdl.get_string_option("ligand");

                parser::FileParser dlpdb (ligand_file,
                                          parser::pdb_read_options::docked_poses_only |
                                          parser::pdb_read_options::skip_atom |
                                          parser::pdb_read_options::all_models);

                dlpdb.parse_molecule(ligand_mols);
                ligand_mols.compute_idatm_type();
                receptor_mols.erase_hydrogen();

                std::map< std::string, std::unique_ptr<score::Score> > scoring_map;

                for (std::string ref : {"complete","reduced"})
                for (std::string comp: {"mean","cumulative"})
                for (std::string func: {"radial","normalized_frequency"})
                for (auto cutoff:{4,5,6,7,8,9,10,11,12,13,14,15}) {
                        std::string score_name = func + "_" + comp + "_" + ref + "_" + to_string(cutoff);
                        scoring_map[score_name] = std::unique_ptr<score::Score>(new score::Score(comp, ref, func,cutoff));

                        scoring_map[score_name]->define_composition(receptor_mols.get_idatm_types(),
                                                 ligand_mols.get_idatm_types())
                        .process_distributions_file(cmdl.get_string_option("dist"))
                        .compile_scoring_function();

                }

                molib::Atom::Grid gridrec(receptor_mols[0].get_atoms());
                std::map< string, vector<double> > output;

                std::vector<std::thread> threads;
                size_t num_threads = cmdl.get_int_option("ncpu");
                std::mutex lock_addition;

                for ( size_t thread_id = 0; thread_id < num_threads; ++thread_id)
                threads.push_back (std::thread ([&, thread_id] {
                        for (size_t i = thread_id; i < ligand_mols.size(); i+=num_threads) {
                        
                                molib::Molecule& ligand  = ligand_mols[i];

                                std::vector<double> results_for_molec;
                                results_for_molec.reserve(96);

                                for (std::string comp: {"mean","cumulative"})
                                for (std::string func: {"radial","normalized_frequency"})
                                for (std::string ref : {"complete","reduced"})
                                for (auto cutoff:{4,5,6,7,8,9,10,11,12,13,14,15}) {
                                        std::string score_name = func + "_" + comp + "_" + ref + "_" + to_string(cutoff);
                                        double new_score = scoring_map[score_name]->non_bonded_energy (gridrec, ligand);
                                        results_for_molec.push_back(new_score);
                                }

                                lock_guard<std::mutex> guard (lock_addition);
                                output.emplace(make_pair(ligand.name(), results_for_molec));
                } } ));

                for (auto &thread : threads) {
                        thread.join();
                }

                receptor_mols.clear();
                ligand_mols.clear();
                for (auto lig : output) {
                        std::cout << lig.first;
                        for (auto score : lig.second)
                                cout << ',' << score;
                        cout << "\n";
                }

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}

