#include "program/target.hpp"

#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "modeler/modeler.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

#include "fileout/fileout.hpp"

using namespace std;
using namespace Program;


// NOTE: Bondorder not included (yet?)
typedef std::map< std::tuple<int,int, size_t >, size_t >         BondCounts;
typedef std::map< std::tuple<int,int,int, size_t >, size_t >     AngleCounts;
typedef std::map< std::tuple<int,int,int,int, size_t >, size_t > DihedralCounts;

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                const double length_bin_size = 0.05;
                const double angle_bin_size  = 0.05;
                const double torsion_bin_size= 0.05;

                for (int i = 1; i < argc; ++i) {
                        Parser::FileParser input(argv[i], Parser::pdb_read_options::all_models);
                        Molib::Molecules mols;
                        input.parse_molecule(mols);
                        
                        mols.compute_idatm_type();
                        mols.erase_hydrogen();
                        
                        auto all_my_atoms = mols[0].get_atoms();
                        Molib::Atom::Graph  graph = Molib::Atom::create_graph(all_my_atoms);
                        auto atom_ring_map = graph.vertex_rings();

                        BondCounts     bc;
                        AngleCounts    ac;
                        DihedralCounts dc;
                        for ( const auto& patom1 : all_my_atoms ) {
                                for ( const auto& atom2 : *patom1 ) {

                                        const Molib::Atom& atom1 = *patom1;

                                        double distance = atom1.crd().distance(atom2.crd());
                                        //size_t bin = std::floor(distance / bin_size);
                                        //auto bond_tuple = std::make_tuple( bond_id.first, bond_id.second, bin);
                                        
                                        /*
                                        if ( !bc.count( bond_tuple ) ) {
                                                bc[bond_tuple] = 1;
                                        } else {
                                                ++bc[bond_tuple];
                                        }*/

                                        cout << "BOND: ";
                                        cout << help::idatm_unmask[atom1.idatm_type()] << " ";
                                        cout << help::idatm_unmask[atom2.idatm_type()] << " ";
                                        cout << distance << endl;

                                        for ( const auto& atom3 : atom2 ) {
                                            
                                                if (atom1.atom_number() == atom3.atom_number())
                                                        continue;
                                            
                                                double angle = Geom3D::angle(atom1.crd(), atom2.crd(), atom3.crd());
                                                /*auto angle_tuple = std::make_tuple( min(atom1->idatm_type(), atom2->idatm_type()),
                                                                                    atom2->idatm_type(),
                                                                                    min(atom1->idatm_type(), atom2->idatm_type()) );*/

                                                cout << "ANGLE: ";
                                                if ( atom1.idatm_type() < atom3.idatm_type() ) {
                                                        cout << help::idatm_unmask[atom1.idatm_type()];
                                                        atom_ring_map[&atom1].size() == 0 ? cout << " " : cout << "_" << atom_ring_map[&atom1].at(0) << " ";
                                                        cout << help::idatm_unmask[atom2.idatm_type()];
                                                        atom_ring_map[&atom2].size() == 0 ? cout << " " : cout << "_" << atom_ring_map[&atom2].at(0) << " ";
                                                        cout << help::idatm_unmask[atom3.idatm_type()];
                                                        atom_ring_map[&atom3].size() == 0 ? cout << " " : cout << "_" << atom_ring_map[&atom3].at(0) << " ";
                                                } else {
                                                        cout << help::idatm_unmask[atom3.idatm_type()];
                                                        atom_ring_map[&atom3].size() == 0 ? cout << " " : cout << "_" << atom_ring_map[&atom3].at(0) << " ";
                                                        cout << help::idatm_unmask[atom2.idatm_type()];
                                                        atom_ring_map[&atom2].size() == 0 ? cout << " " : cout << "_" << atom_ring_map[&atom2].at(0) << " ";
                                                        cout << help::idatm_unmask[atom1.idatm_type()];
                                                        atom_ring_map[&atom1].size() == 0 ? cout << " " : cout << "_" << atom_ring_map[&atom1].at(0) << " ";

                                                }
                                                cout << Geom3D::degrees(angle) << endl;

                                                for ( const auto& atom4 : atom3 ) {
                                                        if (atom4.atom_number() == atom3.atom_number())
                                                                continue;

                                                        double dihedral = Geom3D::dihedral(atom1.crd(), atom2.crd(), atom3.crd(), atom4.crd());

                                                        cout << "DIHEDRAL: ";
                                                        if ( atom1.idatm_type() < atom4.idatm_type() ) {
                                                                cout << help::idatm_unmask[atom1.idatm_type()] << " ";
                                                                cout << help::idatm_unmask[atom2.idatm_type()] << " ";
                                                                cout << help::idatm_unmask[atom3.idatm_type()] << " ";
                                                                cout << help::idatm_unmask[atom4.idatm_type()] << " ";
                                                        } else if ( atom1.idatm_type() > atom4.idatm_type() ) {
                                                                cout << help::idatm_unmask[atom4.idatm_type()] << " ";
                                                                cout << help::idatm_unmask[atom3.idatm_type()] << " ";
                                                                cout << help::idatm_unmask[atom2.idatm_type()] << " ";
                                                                cout << help::idatm_unmask[atom1.idatm_type()] << " ";
                                                        } else {
                                                                if ( atom2.idatm_type() < atom3.idatm_type() ) {
                                                                        cout << help::idatm_unmask[atom1.idatm_type()] << " ";
                                                                        cout << help::idatm_unmask[atom2.idatm_type()] << " ";
                                                                        cout << help::idatm_unmask[atom3.idatm_type()] << " ";
                                                                        cout << help::idatm_unmask[atom4.idatm_type()] << " ";
                                                                } else { // does not matter if they are equal, it will be the same result
                                                                        cout << help::idatm_unmask[atom4.idatm_type()] << " ";
                                                                        cout << help::idatm_unmask[atom3.idatm_type()] << " ";
                                                                        cout << help::idatm_unmask[atom2.idatm_type()] << " ";
                                                                        cout << help::idatm_unmask[atom1.idatm_type()] << " ";
                                                                }
                                                        }
                                                        cout << dihedral << "\n";
                                                }
                                        }
                                }
                        }
                }


        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}


