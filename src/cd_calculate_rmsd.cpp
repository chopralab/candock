#include "program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Benchmark main_timer;
                main_timer.display_time("Starting");

                if ( argc <= 1 ) {
                        cerr << "You MUST supply an argument!" << endl;
                        return 1;
                }

                Parser::FileParser mol1(argv[1], Parser::pdb_read_options::docked_poses_only | Parser::pdb_read_options::skip_atom | Parser::pdb_read_options::all_models);

                Molib::Molecules mols1;

                mol1.parse_molecule(mols1);

                std::ofstream output_matrix("outmatrix.tbl");

                if (mols1.size() <= 0) {
                        output_matrix.close();
                        return 0;
                }

                if (mols1.size() == 1) {
                        output_matrix << 0 << endl;
                        output_matrix.close();
                        return 0;
                }

                std::vector< std::vector<double> > previous(mols1.size());

                for (size_t i = 0; i < mols1.size(); ++i) {

                        for (size_t l = 0; l < i; ++l)
                                output_matrix << setw (12) << setprecision (8) << previous[l].at (i - (1 + l));

                        output_matrix << setw (12) << setprecision (8) << 0;

                        for (size_t j = i + 1; j < mols1.size(); ++j) {
                                double rmsd = mols1[i].compute_rmsd_ord(mols1[j]);
                                previous[i].push_back (rmsd);
                                output_matrix << setw (12) << setprecision (8) << rmsd;
                        }

                        output_matrix << '\n';
                }

                output_matrix.close();

                main_timer.display_time("Finished");

        } catch (exception& e) {
                log_error << e.what() << endl;
                return 1;
        }
        return 0;
}
