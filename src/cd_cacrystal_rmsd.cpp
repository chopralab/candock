#include "program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Benchmark main_timer;
                main_timer.display_time("Starting");

                if ( argc <= 2 ) {
                        cerr << "You MUST supply at least two arguments!" << endl;
                        return 1;
                }

                Parser::FileParser crys(argv[1], Parser::pdb_read_options::first_model);
                Parser::FileParser mol1(argv[2], Parser::pdb_read_options::docked_poses_only | Parser::pdb_read_options::skip_atom | Parser::pdb_read_options::all_models);

                Molib::Molecules cryst;
                Molib::Molecules mols1;

                crys.parse_molecule(cryst);
                mol1.parse_molecule(mols1);

                std::ofstream output_matrix;
                
                if (argc >= 4) {
                        output_matrix.open(argv[3]);
                } else {
                        output_matrix.open("rmsds.lst");
                }

                if (mols1.size() <= 0) {
                        output_matrix.close();
                        return 0;
                }

                cryst[0].erase_properties();
                Molib::Atom::Graph cryst_graph =  Molib::Atom::create_graph(cryst[0].get_atoms());
                std::vector< Molib::Atom::Graph > atom_graphs;

                for ( Molib::Molecule& mol : mols1 ) {
                        Molib::Atom::Vec atoms = mol.get_atoms();
                        
                        int reenum = 0;
                        for ( auto &atom : atoms ) {
                                atom->set_atom_number(++reenum);
                                atom->erase_properties();
                        }

                        atom_graphs.push_back( Molib::Atom::create_graph(atoms));
                }

                for (size_t i = 0; i < atom_graphs.size(); ++i) {

                        double rmsd = Molib::Atom::compute_rmsd(cryst_graph, atom_graphs[i]);
                        output_matrix << rmsd << '\n';
                }

                output_matrix.close();

                main_timer.display_time("Finished");

        } catch (exception& e) {
                log_error << e.what() << endl;
                return 1;
        }
        return 0;
}

