#include "program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;

int main(int argc, char* argv[]) {
        try {

                Benchmark main_timer;
                main_timer.display_time("Starting");

                Parser::FileParser mol1(argv[1], Parser::pdb_read_options::all_models | Parser::pdb_read_options::skip_atom);

                Molib::Molecules mols1;

                mol1.parse_molecule(mols1);

                std::vector< std::vector<double> > previous(mols1.size());

                std::ofstream output_matrix("outmaxtrix.tbl");

                if (mols1.size() <= 1) {
                        output_matrix << 0 << endl;
                } else {
                        for (size_t i = 0; i < mols1.size(); ++i) {

                                for (size_t l = 0; l < i; ++l)
                                        output_matrix << setw (12) << setprecision (8) << previous[l].at (i - (1 + l));

                                output_matrix << setw (12) << setprecision (8) << 0;

                                for (size_t j = i + 1; j < mols1.size(); ++j) {
                                        double rmsd = mols1[i].compute_rmsd (mols1[j]);
                                        previous[i].push_back (rmsd);
                                        output_matrix << setw (12) << setprecision (8) << rmsd;
                                }

                                output_matrix<< endl;
                        }                
                }

                output_matrix.close();

                main_timer.display_time("Finished");

        } catch (exception& e) {
                log_error << e.what() << endl;
                return 1;
        }
        return 0;
}
