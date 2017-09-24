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


class MoleculeBondExtractor {

public:
        typedef std::tuple< const Molib::Atom*, const Molib::Atom* > BondStretch;
        typedef std::map< BondStretch, double> BondStretches;
        typedef std::tuple< const Molib::Atom*, const Molib::Atom*,
                            const Molib::Atom* > BondAngle;
        typedef std::map< BondAngle, double> BondAngles;
        typedef std::tuple< const Molib::Atom*, const Molib::Atom*,
                            const Molib::Atom*, const Molib::Atom* > BondDihedral;
        typedef std::map< BondDihedral, double> BondDihedrals;

private:
        BondStretches bs;
        BondAngles    ba;
        BondDihedrals bd;
        Glib::Graph<Molib::Atom>::VertexRingMap all_rings;
        std::map<const Molib::Atom*, size_t> substitutions;

        std::string __print_ring_info( const std::vector<size_t>& ring_vec) const {
                if (ring_vec.size() == 0) {
                        return "0";
                }
                if (ring_vec.size() > 1) {
                        return "BR";
                }
                return std::to_string(ring_vec.at(0));
        }
        
public:

        void addStretch ( BondStretch stretch ) {
                double distance = get<0>(stretch)->crd().distance(get<1>(stretch)->crd());
                bs[stretch] = distance;
        }

        void addAngle ( BondAngle angle) {
                const Molib::Atom *atom1, *atom2, *atom3;
                std::tie(atom1, atom2, atom3) = angle;

                if (atom1->atom_number() == atom3->atom_number())
                        return;

                double angle_val = Geom3D::angle(atom1->crd(), atom2->crd(), atom3->crd());
                if ( atom1->idatm_type() < atom3->idatm_type() ) {
                        ba[make_tuple(atom1, atom2, atom3)] = angle_val;
                } else {
                        ba[make_tuple(atom3, atom2, atom1)] = angle_val;
                }
        }

        void addDihedral (BondDihedral dihedral_atoms) {

                const Molib::Atom *atom1, *atom2, *atom3, *atom4;
                std::tie(atom1, atom2, atom3, atom4) = dihedral_atoms;
            
                if (atom4->atom_number() == atom3->atom_number())
                        return;

                double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                                   atom3->crd(), atom4->crd());

                if ( atom1->idatm_type() < atom4->idatm_type() ) {
                        bd[make_tuple(atom1, atom2, atom3, atom4)] = dihedral;
                } else if ( atom1->idatm_type() > atom4->idatm_type() ) {
                        bd[make_tuple(atom4, atom3, atom2, atom1)] = dihedral;
                } else {
                        if ( atom2->idatm_type() < atom3->idatm_type() ) {
                                bd[make_tuple(atom1, atom2, atom3, atom4)] = dihedral;
                        } else { // does not matter if they are equal, it will be the same result
                                bd[make_tuple(atom4, atom3, atom2, atom1)] = dihedral;
                        }
                }
        }

        void addMolecule( const Molib::Molecule& mol ) {
                auto all_my_atoms = mol.get_atoms();
                Molib::Atom::Graph  graph = Molib::Atom::create_graph(all_my_atoms);
                auto atom_ring_map = graph.vertex_rings();
                all_rings.insert( atom_ring_map.begin(), atom_ring_map.end());

                for ( const auto& patom1 : all_my_atoms ) {
                        for ( const auto& atom2 : *patom1 ) {

                                const Molib::Atom& atom1 = *patom1;

                                addStretch(make_tuple(&atom1, &atom2));

                                for ( const auto& atom3 : atom2 ) {
                                        addAngle(make_tuple(&atom1, &atom2, &atom3));
                                        for ( const auto& atom4 : atom3 ) {
                                                addDihedral(make_tuple(&atom1, &atom2,
                                                                       &atom3, &atom4));
                                        }
                                }
                        }

                        substitutions[patom1] = patom1->size() - patom1->get_num_hydrogens();
                }
        }

        void printBonds () const {
                for ( const auto& bond : bs ) {
                        cout << "BOND: ";
                        cout << help::idatm_unmask[ get<0>(bond.first)->idatm_type() ] << " ";
                        cout << help::idatm_unmask[ get<1>(bond.first)->idatm_type() ] << " ";
                        cout << bond.second;
                        cout << "\n";
                }
        }

        void printAngles () {
                for ( const auto& angle : ba ) {
                        cout << "ANGLE: ";
                        cout << help::idatm_unmask[ get<0>(angle.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<0>(angle.first)]) << " ";
                        cout << help::idatm_unmask[ get<1>(angle.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<1>(angle.first)]) << " ";
                        cout << help::idatm_unmask[ get<2>(angle.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<2>(angle.first)]) << " ";
                        cout << Geom3D::degrees(angle.second);
                        cout << "\n";
                }
        }
        
        void printDihedrals () {
                for ( const auto& dihedral : bd ) {
                        cout << "DIHEDRAL: ";
                        cout << help::idatm_unmask[ get<0>(dihedral.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<0>(dihedral.first)]) << "_";
                        cout << substitutions[get<0>(dihedral.first)] << " ";
                        cout << help::idatm_unmask[ get<1>(dihedral.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<1>(dihedral.first)]) << "_";
                        cout << substitutions[get<1>(dihedral.first)] << " ";
                        cout << help::idatm_unmask[ get<2>(dihedral.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<2>(dihedral.first)]) << "_";
                        cout << substitutions[get<2>(dihedral.first)] << " ";
                        cout << help::idatm_unmask[ get<3>(dihedral.first)->idatm_type() ] << "_";
                        cout << __print_ring_info(all_rings[get<3>(dihedral.first)]) << "_";
                        cout << substitutions[get<3>(dihedral.first)] << " ";
                        cout << Geom3D::degrees(dihedral.second);
                        cout << "\n";
                }
        }


};

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                MoleculeBondExtractor MBE;

                for (int i = 1; i < argc; ++i) {
                        Parser::FileParser input(argv[i], Parser::pdb_read_options::all_models);
                        Molib::Molecules mols;
                        input.parse_molecule(mols);
                        
                        mols.compute_idatm_type();
                        mols.erase_hydrogen();
                        
                        MBE.addMolecule(mols[0]);
                }

                MBE.printBonds();
                MBE.printAngles();
                MBE.printDihedrals();

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}
