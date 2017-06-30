#include "modeler/forcefield.hpp"
#include "modeler/modeler.hpp"
#include "program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {
#ifndef _WINDOWS

                Inout::Logger::set_all_stderr(true);

                help::Options::set_options( new Program::CmdLnOpts( 
                    argc, argv));

                if ( argc <= 1 ) {
                        cerr << "You supply an argument!" << endl;
                        return 1;
                }

                Parser::FileParser mol1(argv[1], Parser::pdb_read_options::docked_poses_only | Parser::pdb_read_options::skip_atom | Parser::pdb_read_options::all_models);

                Molib::Molecules mols1;
                mol1.parse_molecule(mols1);

                if (mols1.size() <= 0) {
                        return 0;
                }

                size_t num_threads = 1;
                if (argc >= 3) {
                        num_threads = std::atoi(argv[2]);
                }

                OMMIface::ForceField ffield;

                ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
                      .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                      .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                ffield.insert_topology (mols1[0]);

                std::vector<double> potentials(mols1.size());

                std::vector<std::thread> threads;

                for ( size_t thread_id = 0; thread_id < num_threads; ++thread_id)
                threads.push_back( (std::thread ([&, thread_id] {
                        for ( size_t i = thread_id; i < mols1.size(); i+=num_threads) {

                                OMMIface::Modeler modeler (ffield,
                                                           cmdl.get_string_option("fftype"), cmdl.get_int_option("cutoff"),
                                                           cmdl.get_double_option("mini_tol"), cmdl.get_int_option("max_iter"),
                                                           cmdl.get_int_option("update_freq"), cmdl.get_double_option("pos_tol"),
                                                           false, cmdl.get_double_option("dynamic_step_size"),
                                                           cmdl.get_double_option("temperature"), cmdl.get_double_option("friction"));

                                modeler.add_topology (mols1[i].get_atoms());

                                modeler.init_openmm();
                                modeler.add_crds (mols1[i].get_atoms(), mols1[i].get_crds());
                                modeler.init_openmm_positions();

                                const double potential_energy = modeler.potential_energy();

                                potentials[i] = potential_energy;
                        }
                })));
                
                for (auto &thread : threads) {
                        thread.join();
                }

                for ( auto &poten : potentials ) {
                        cout << poten << endl;
                }
#endif
        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}


