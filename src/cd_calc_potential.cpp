#include "modeler/forcefield.hpp"
#include "modeler/modeler.hpp"
#include "program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                help::Options::set_options( new Program::CmdLnOpts( 
                    argc, argv));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                if ( argc <= 1 ) {
                        cerr << "You supply an argument!" << endl;
                        return 1;
                }

                Parser::FileParser mol1(argv[1], Parser::pdb_read_options::docked_poses_only | Parser::pdb_read_options::skip_atom | Parser::pdb_read_options::all_models);

                Molib::Molecules mols1;

                mol1.parse_molecule(mols1);

                std::ofstream output_matrix;
                
                if (argc >= 3) {
                        output_matrix.open(argv[2]);
                } else {
                        output_matrix.open("poten.lst");
                }
                
                if (mols1.size() <= 0) {
                        output_matrix.close();
                        return 0;
                }

                OMMIface::ForceField ffield;

                ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
                      .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                      .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                ffield.insert_topology (mols1[0]);

                for ( size_t i = 0; i < mols1.size(); ++i) {

                        OMMIface::Modeler modeler (ffield, cmdl.get_string_option("fftype"), cmdl.get_int_option("cutoff"),
                                           cmdl.get_double_option("mini_tol"), cmdl.get_int_option("max_iter"), cmdl.get_int_option("update_freq"), 
                                           cmdl.get_double_option("pos_tol"), false, 2.0);

                        modeler.add_topology (mols1[i].get_atoms());

                        modeler.init_openmm();
                        modeler.add_crds (mols1[i].get_atoms(), mols1[i].get_crds());
                        modeler.init_openmm_positions();

                        const double potential_energy = modeler.potential_energy();

                        output_matrix << potential_energy << '\n';
                }
                
                output_matrix.close();

                main_timer.display_time("Finished");

        } catch (exception& e) {
                log_error << e.what() << endl;
                return 1;
        }
        return 0;
}


