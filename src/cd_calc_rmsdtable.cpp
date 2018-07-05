#include "candock/program/fragmentligands.hpp"
#include "version.hpp"

////////////////// FRAGMENTING OF LIGANDS ///////////////////////////

using namespace std;
using namespace candock;

int main(int argc, char* argv[]) {
        try {
#ifndef _WINDOWS
                if ( argc <= 1 ) {
                        cerr << "You MUST supply an argument!" << endl;
                        cerr << "A fourth argument gives the number of threads (defaults to 1)" << endl;
                        return 1;
                }

                Inout::Logger::set_all_stderr(true);

                parser::FileParser mol1(argv[1], parser::pdb_read_options::docked_poses_only |
                                                 parser::pdb_read_options::skip_atom |
                                                 parser::pdb_read_options::all_models);

                molib::Molecules mols1;

                mol1.parse_molecule(mols1);

                if (mols1.size() <= 0) {
                        return 0;
                }

                if (mols1.size() == 1) {
                        std::cout << 0 << endl;
                        return 0;
                }


                size_t num_threads = 1;
                if ( argc >= 3 ) {
                        num_threads = std::atoi(argv[2]);
                }

                std::vector< std::vector<double> > matrix(mols1.size(), std::vector<double>(mols1.size()));

                std::vector<std::thread> threads;

                for ( size_t thread_id = 0; thread_id < num_threads; ++thread_id)
                threads.push_back (std::thread ([&, thread_id] {
                        for (size_t i = thread_id; i < mols1.size(); i += num_threads) {
                                for (size_t j = i + 1; j < mols1.size(); ++j) {
                                        double rmsd = mols1[i].compute_rmsd_ord(mols1[j]);
                                        matrix[i][j] = rmsd;
                                        matrix[j][i] = rmsd;
                                }

                        }
                }));

                for (auto &thread : threads) {
                        thread.join();
                }


                for (const auto &row : matrix) {
                        for (const auto &col : row) {
                                std::cout << setw (12) << setprecision (8) << col;
                        }
                        cout << '\n';
                }
#endif
        } catch (exception& e) {
                cerr << e.what() << endl;
                return 1;
        }
        return 0;
}
