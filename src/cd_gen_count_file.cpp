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
        typedef std::map< BondDihedral, double> BondImproper;

private:
        BondStretches bs;
        BondAngles    ba;
        BondDihedrals bd;
        BondDihedrals bi; // Impropers

        Glib::Graph<Molib::Atom>::VertexRingMap all_rings;
        std::map<const Molib::Atom*, size_t> substitutions;

        int __comp_atoms ( const Molib::Atom* atom1, const Molib::Atom* atom2) const {

                if ( atom1 == atom2 )
                        return 0;

                if ( atom1->idatm_type() < atom2->idatm_type() ) {
                        return -2;
                } else if ( atom1->idatm_type() > atom2->idatm_type() ) {
                        return 2;
                }
                
                //IDATM types are equal, compare addresses
                if ( atom1 < atom2 ) {
                        return -1;
                } else if ( atom1 > atom2 ) {
                        return 1;
                }
                
                return 0;
        }

        std::string __print_ring_info( const std::vector<size_t>& ring_vec) const {
                if (ring_vec.size() == 0) {
                        return "0";
                }
                if (ring_vec.size() > 1) {
                        return "BR";
                }
                return std::to_string(ring_vec.at(0));
        }

        void __print_atom( const Molib::Atom* atom, std::ostream& os, char delim ) const {
                os << help::idatm_unmask[ atom->idatm_type() ] << delim;
                os << __print_ring_info(all_rings.at(atom)) << delim;
                os << substitutions.at(atom) << delim;

        }

public:

        void addStretch ( BondStretch stretch_atoms ) {

                const Molib::Atom *atom1, *atom2;
                std::tie(atom1, atom2) = stretch_atoms;

                BondStretch s_new;

                if ( atom1 < atom2 ) {
                        s_new = make_tuple(atom1, atom2);
                } else {
                        s_new = make_tuple(atom2, atom1);
                }

                if ( bs.count(s_new) )
                        return;

                double distance = atom1->crd().distance(atom2->crd());
                bs[s_new] = distance;
        }

        void addAngle ( BondAngle angle_atoms) {

                const Molib::Atom *atom1, *atom2, *atom3;
                std::tie(atom1, atom2, atom3) = angle_atoms;

                if (atom1 == atom3)
                        return;

                BondAngle a_new;
                if ( atom1 < atom3 ) {
                        a_new = make_tuple(atom1, atom2, atom3);
                } else {
                        a_new = make_tuple(atom3, atom2, atom1);
                }

                if ( ba.count(a_new) )
                        return;

                double angle = Geom3D::angle(atom1->crd(), atom2->crd(), atom3->crd());
                ba[a_new] = angle;
        }

        void addDihedral (BondDihedral dihedral_atoms) {

                const Molib::Atom *atom1, *atom2, *atom3, *atom4;
                std::tie(atom1, atom2, atom3, atom4) = dihedral_atoms;
            
                if (atom4 == atom2 || atom4 == atom1 || atom3 == atom1)
                        return;


                BondDihedral d_new;

                if ( atom1 < atom4 ) {
                        d_new = make_tuple(atom1, atom2, atom3, atom4);
                } else if ( atom4 < atom1 ) {
                        d_new = make_tuple(atom4, atom3, atom2, atom1);
                }

                if ( bd.count(d_new) )
                        return;

                double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                                   atom3->crd(), atom4->crd());
                bd[d_new] = dihedral;
        }
        
        void addImproper (BondDihedral improper_atoms) {

                const Molib::Atom *atom1, *atom2, *atom3, *atom4;
                std::tie(atom1, atom2, atom3, atom4) = improper_atoms;
            
                if (atom4 == atom3 || atom4 == atom1 || atom3 == atom1)
                        return;

                // Order should be:
                // Lowest, 2nd Lowest, Atom2, highest
                BondDihedral i_new;
                
                const Molib::Atom* highest = max( max( atom1, atom3), atom4);
                const Molib::Atom* lowest  = min( min( atom1, atom3), atom4);
                const Molib::Atom* middle  = max( min(atom1,atom3), min(max(atom1,atom3),atom4));

                i_new = make_tuple (lowest, middle, atom2, highest);

                if ( bd.count(i_new) )
                        return;

                double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                                   atom3->crd(), atom4->crd());
                bd[i_new] = dihedral;
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
                                        for ( const auto& atom4 : atom2 ) {
                                                addImproper(make_tuple(&atom1, &atom2,
                                                                       &atom3, &atom4));
                                        }
                                }
                        }

                        substitutions[patom1] = patom1->size() - patom1->get_num_hydrogens();
                }
        }

        void printBonds (std::ostream& os) const {
                for ( const auto& bond : bs ) {
                        __print_atom(get<0>(bond.first), os, ',');
                        __print_atom(get<1>(bond.first), os, ',');
                        os << bond.second;
                        os << "\n";
                }
        }

        void printAngles (std::ostream& os) const {
                for ( const auto& angle : ba ) {
                        __print_atom(get<0>(angle.first), os, ',');
                        __print_atom(get<1>(angle.first), os, ',');
                        __print_atom(get<2>(angle.first), os, ',');
                        os << Geom3D::degrees(angle.second);
                        os << "\n";
                }
        }
        
        void printDihedrals (std::ostream& os) const {
                for ( const auto& dihedral : bd ) {
                        __print_atom(get<0>(dihedral.first), os, ',');
                        __print_atom(get<1>(dihedral.first), os, ',');
                        __print_atom(get<2>(dihedral.first), os, ',');
                        __print_atom(get<3>(dihedral.first), os, ',');
                        os << Geom3D::degrees(dihedral.second);
                        os << "\n";
                }
        }

        void printImpropers (std::ostream& os) const {
                for ( const auto& improper : bi ) {
                        __print_atom(get<0>(improper.first), os, ',');
                        __print_atom(get<1>(improper.first), os, ',');
                        __print_atom(get<2>(improper.first), os, ',');
                        __print_atom(get<3>(improper.first), os, ',');
                        os << Geom3D::degrees(improper.second);
                        os << "\n";
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

                MBE.printBonds(std::cout);
                MBE.printAngles(std::cout);
                MBE.printDihedrals(std::cout);
                MBE.printImpropers(std::cout);

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}