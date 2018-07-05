#include "candock/program/cmdlnopts.hpp"
#include "candock/program/fragmentligands.hpp"
#include "candock/program/target.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"

#include "candock/score/xscore.hpp"
#include "candock/modeler/forcefield.hpp"

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

                Parser::FileParser drpdb (cmdl.get_string_option("receptor"),
                                          Parser::pdb_read_options::protein_poses_only |
                                          Parser::pdb_read_options::all_models);

                drpdb.parse_molecule(receptor_mols);

                // Check to see if the user set the ligand option, use the receptor as a last resort.
                const std::string &ligand_file = Inout::file_size(cmdl.get_string_option("ligand")) == 0 ?
                                                    cmdl.get_string_option("receptor") :
                                                    cmdl.get_string_option("ligand");

                Parser::FileParser dlpdb (ligand_file,
                                          Parser::pdb_read_options::docked_poses_only |
                                          Parser::pdb_read_options::skip_atom |
                                          Parser::pdb_read_options::all_models);

                dlpdb.parse_molecule(ligand_mols);

                if (receptor_mols.size() == 0 || ligand_mols.size() == 0) {
                        return 0;
                }

                std::vector<std::thread> threads;

                size_t num_threads = cmdl.get_int_option("ncpu");

                std::vector<std::array<double,5>> output;

                output.reserve(ligand_mols.size());

                for ( size_t thread_id = 0; thread_id < num_threads; ++thread_id)
                threads.push_back (std::thread ([&, thread_id] {
                        for (size_t i = thread_id; i < receptor_mols.size(); i+=num_threads) {
                        
                                molib::Molecule& protein = receptor_mols[i];
                                molib::Molecule& ligand  = ligand_mols[i];

                                molib::Atom::Grid gridrec(protein.get_atoms());                        

                                output.emplace_back(score::vina_xscore(gridrec, ligand.get_atoms()));
                } } ));

                for (auto &thread : threads) {
                        thread.join();
                }

                receptor_mols.clear();
                ligand_mols.clear();
                for (const auto& values : output ) {
                        for (const auto& component : values) {
                                std::cout << component << ",";
                        }
                        std::cout << "\n";
                }

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}

