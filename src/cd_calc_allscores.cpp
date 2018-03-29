#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

#include "score/score.hpp"
#include "modeler/forcefield.hpp"

#include <memory>

using namespace std;
using namespace Program;

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

                Molib::Molecules receptor_mols;
                Molib::Molecules ligand_mols;

                Parser::FileParser drpdb (cmdl.get_string_option("receptor"),
                                          Parser::pdb_read_options::protein_poses_only |
                                          Parser::pdb_read_options::all_models);

                drpdb.parse_molecule(receptor_mols);
                receptor_mols.compute_idatm_type();
                receptor_mols.erase_hydrogen();

                // Check to see if the user set the ligand option, use the receptor as a last resort.
                const std::string &ligand_file = Inout::file_size(cmdl.get_string_option("ligand")) == 0 ?
                                                    cmdl.get_string_option("receptor") :
                                                    cmdl.get_string_option("ligand");

                Parser::FileParser dlpdb (ligand_file,
                                          Parser::pdb_read_options::docked_poses_only |
                                          Parser::pdb_read_options::skip_atom |
                                          Parser::pdb_read_options::all_models);

                dlpdb.parse_molecule(ligand_mols);
                ligand_mols.compute_idatm_type();
                receptor_mols.erase_hydrogen();

                std::map< std::string, std::unique_ptr<Score::Score> > scoring_map;

                for (std::string ref : {"complete","reduced"})
                for (std::string comp: {"mean","cumulative"})
                for (std::string func: {"radial","normalized_frequency"})
                for (auto cutoff:{4,5,6,7,8,9,10,11,12,13,14,15}) {
                        std::string score_name = func + "_" + comp + "_" + ref + "_" + to_string(cutoff);
                        scoring_map[score_name] = std::make_unique<Score::Score>(ref, comp,func,cutoff, cmdl.get_double_option("step"));

                        scoring_map[score_name]->define_composition(receptor_mols.get_idatm_types(),
                                                 ligand_mols.get_idatm_types())
                        .process_distributions_file(cmdl.get_string_option("dist"))
                        .compile_scoring_function();

                }

                Molib::Atom::Grid gridrec(receptor_mols[0].get_atoms());
                std::map< string, vector<double> > output;

                std::vector<std::thread> threads;
                size_t num_threads = cmdl.get_int_option("ncpu");
                std::mutex lock_addition;

                for ( size_t thread_id = 0; thread_id < num_threads; ++thread_id)
                threads.push_back (std::thread ([&, thread_id] {
                        for (size_t i = thread_id; i < ligand_mols.size(); i+=num_threads) {
                        
                                Molib::Molecule& ligand  = ligand_mols[i];

                                std::vector<double> results_for_molec;
                                results_for_molec.reserve(96);

                                for (std::string comp: {"mean","cumulative"})
                                for (std::string func: {"radial","normalized_frequency"})
                                for (std::string ref : {"complete","reduced"})
                                for (auto cutoff:{4,5,6,7,8,9,10,11,12,13,14,15}) {
                                        std::string score_name = func + "_" + comp + "_" + ref + "_" + to_string(cutoff);
                                        results_for_molec.push_back(scoring_map[score_name]->non_bonded_energy (gridrec, ligand));
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

