#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

#include "score/score.hpp"
#include "modeler/forcefield.hpp"

using namespace std;
using namespace Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

void score_thread(Parser::FileParser &drpdb, Parser::FileParser &dlpdb, const Score::Score& score,
        std::vector<double> &output, std::mutex& mutex, int start_index) {


}

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

                // Check to see if the user set the ligand option, use the receptor as a last resort.
                const std::string &ligand_file = Inout::file_size(cmdl.get_string_option("ligand")) == 0 ?
                                                    cmdl.get_string_option("receptor") :
                                                    cmdl.get_string_option("ligand");

                Parser::FileParser dlpdb (ligand_file,
                                          Parser::pdb_read_options::docked_poses_only |
                                          Parser::pdb_read_options::skip_atom |
                                          Parser::pdb_read_options::all_models);

                dlpdb.parse_molecule(ligand_mols);

                Score::Score score(cmdl.get_string_option("ref"), cmdl.get_string_option("comp"),
                                   cmdl.get_string_option("func"),cmdl.get_int_option("cutoff"),
                                   cmdl.get_double_option("step"));

                score.define_composition(receptor_mols.get_idatm_types(),
                                         ligand_mols.get_idatm_types())
                     .process_distributions_file(cmdl.get_string_option("dist"))
                     .compile_scoring_function();

                std::vector<std::thread> threads;

                size_t num_threads = cmdl.get_int_option("ncpu");
                std::vector<double> output (receptor_mols.size());
                for ( size_t thread_id = 0; thread_id < num_threads; ++thread_id)
                threads.push_back (std::thread ([&, thread_id] {
                        for (size_t i = thread_id; i < receptor_mols.size(); i+=num_threads) {
                        
                                Molib::Molecule& protein = receptor_mols[i];
                                Molib::Molecule& ligand  = ligand_mols[i];

                                Molib::Atom::Grid gridrec(protein.get_atoms());                        

                                output[i] = score.non_bonded_energy (gridrec, ligand);
                } } ));

                for (auto &thread : threads) {
                        thread.join();
                }

                receptor_mols.clear();
                ligand_mols.clear();
                for (size_t i = 0; i < output.size(); ++i) {
                        std::cout << output[i] << '\n';
                }

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}

