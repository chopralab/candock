#include "program/target.hpp"

#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "modeler/modeler.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

#include "fileout/fileout.hpp"

using namespace std;
using namespace Program;

// TODO: Break this into two classes, one to extract data from a molecule, another to bin
class MoleculeBondExtractor {

public:

        // typedefs for raw Bond information
        typedef std::tuple< const Molib::Atom*,
                            const Molib::Atom* > BondStretch;
        typedef std::map< BondStretch, double>   BondStretches;
        typedef std::tuple< const Molib::Atom*,
                            const Molib::Atom*,
                            const Molib::Atom* > BondAngle;
        typedef std::map< BondAngle, double>     BondAngles;
        typedef std::tuple< const Molib::Atom*,
                            const Molib::Atom*,
                            const Molib::Atom*,
                            const Molib::Atom* > BondDihedral;
        typedef std::map< BondDihedral, double>  BondDihedrals;

        // typedefs for binned Bond information
        // NOTE: Bondorder not included (yet?)
        typedef std::tuple<int,int, size_t>                          StretchBin;
        typedef std::map< StretchBin, size_t >                       StretchCounts;
        typedef std::tuple<int,int,int,int,int,int, size_t >         AngleBin;
        typedef std::map< AngleBin, size_t >                         AngleCounts;
        typedef std::tuple<int,int,int,int,
                           int,int,int,int,
                           int,int,int,int,int >                  DihedralBin;
        typedef std::map< DihedralBin, size_t >                      DihedralCounts;

private:
        BondStretches bs;
        BondAngles    ba;
        BondDihedrals bd;
        BondDihedrals bi; // Impropers

        StretchCounts  bsc;
        AngleCounts    bac;
        DihedralCounts bdc;
        DihedralCounts bic;
        
        Glib::Graph<Molib::Atom>::Cycles all_rings;
        Glib::Graph<Molib::Atom>::VertexRingMap all_rings_sizes;
        std::map<const Molib::Atom*, size_t> substitutions;

        int __comp_atoms ( const Molib::Atom* atom1, const Molib::Atom* atom2) const {

                if ( atom1 == atom2 )
                        return 0;

                if ( atom1->idatm_type() < atom2->idatm_type() ) {
                        return -4;
                } else if ( atom1->idatm_type() > atom2->idatm_type() ) {
                        return 4;
                }

                if ( __ring_size(atom1) < __ring_size(atom2) ) {
                        return -3;
                } else if (__ring_size(atom1) > __ring_size(atom2)) {
                        return 3;
                }

                if ( substitutions.at(atom1) < substitutions.at(atom2) ) {
                        return -2;
                } else if (substitutions.at(atom1) > substitutions.at(atom2)) {
                        return 2;
                }
                                
                return 0;
        }

        int __ring_size( const Molib::Atom* atom ) const {
                const std::vector<size_t>& ring_vec = all_rings_sizes.at(atom);
                if (ring_vec.size() == 0) { // Not in a ring
                        return 0;
                }
                if (ring_vec.size() > 1) { // Bridged or Spiro

                        if ( atom->element() != Molib::Element::C) // Spiro must be carbon
                                return 2;
                        
                        if ( substitutions.at(atom) != 4 ) { // Spiro must have four non-hydrogen neighbors
                                return 2;
                        }

                        for ( const Molib::Bond* bond : atom->get_bonds() ) {
                                const Molib::Atom& other = bond->second_atom(*atom);

                                if ( all_rings_sizes.at(&other).size() == 0 ) { // Spiro must be connected to only rings
                                        return 2;
                                }
                                
                                if ( all_rings_sizes.at(&other).size() == 1 ) { // Vast majority of cases present here
                                        continue;
                                }

                                // Missing quite a few Spiro cases here, they can be added later.
                                size_t shared_rings = 0;
                                for ( const auto& thing : all_rings ) {
                                        // Are they in the same ring?
                                        if ( thing.count( const_cast<Molib::Atom*>(atom) ) &&
                                             thing.count( const_cast<Molib::Atom*>(&other) ) ) {
                                                shared_rings++;
                                        }
                                }

                                // FIXME: should be adjusted for large ring structures.
                                if ( shared_rings >= 2 )
                                        return 2;
                        }
                        return 1;
                }
                return ring_vec.at(0);
        }

        void __print_atom( const Molib::Atom* atom, std::ostream& os, char delim ) const {
                os << help::idatm_unmask[ atom->idatm_type() ] << delim;
                os << __ring_size(atom) << delim;
                os << substitutions.at(atom) << delim;

        }

public:

        void reset() {
                substitutions.clear();
                all_rings.clear();
                all_rings_sizes.clear();
        }

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
                const Molib::Atom* middle  = max( min( atom1, atom3), min(max(atom1,atom3),atom4));

                i_new = make_tuple (lowest, middle, atom2, highest);

                if ( bi.count(i_new) )
                        return;

                double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                                   atom3->crd(), atom4->crd());
                bi[i_new] = dihedral;
        }


        bool addMolecule( const Molib::Molecule& mol ) {

                auto all_my_atoms = mol.get_atoms();
                
                if (!all_my_atoms.size()) {
                        return false;
                }
                
                Molib::Atom::Graph  graph = Molib::Atom::create_graph(all_my_atoms);
                auto atom_ring_map = graph.vertex_rings();
                all_rings_sizes.insert( atom_ring_map.begin(), atom_ring_map.end());
                auto atom_all_rings= graph.find_rings();
                all_rings.insert(atom_all_rings.begin(), atom_all_rings.end());

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
                
                return true;
        }

        void binStretches( const double bond_bin_size ) {
                 for ( const auto& bond : bs ) {
                        const size_t bin = floor(bond.second / bond_bin_size);
                        
                        StretchBin sb;
                        if ( get<0>(bond.first)->idatm_type() <
                             get<1>(bond.first)->idatm_type() ) {
                                sb = make_tuple( get<0>(bond.first)->idatm_type(),
                                                 get<1>(bond.first)->idatm_type(),
                                                 bin
                                               );
                        } else {
                                sb = make_tuple( get<1>(bond.first)->idatm_type(),
                                                 get<0>(bond.first)->idatm_type(),
                                                 bin
                                               );
                        }
                        
                        if ( ! bsc.count(sb) )
                            bsc[sb] = 1;
                        else
                            ++bsc[sb];
                 }
                 bs.clear();
        }

        void binAngles( const double angle_bin_size ) {
                 for ( const auto& angle : ba ) {
                        const size_t bin = floor(angle.second / angle_bin_size);

                        AngleBin ab;
                        if ( __comp_atoms(get<0>(angle.first), get<2>(angle.first) ) < 0 ) {
                                ab = make_tuple( get<0>(angle.first)->idatm_type(),
                                                 __ring_size(get<0>(angle.first)),
                                                 get<1>(angle.first)->idatm_type(),
                                                 __ring_size(get<1>(angle.first)),
                                                 get<2>(angle.first)->idatm_type(),
                                                 __ring_size(get<2>(angle.first)),
                                                 bin
                                               );
                        } else {
                                ab = make_tuple( get<2>(angle.first)->idatm_type(),
                                                 __ring_size(get<2>(angle.first)),
                                                 get<1>(angle.first)->idatm_type(),
                                                 __ring_size(get<1>(angle.first)),
                                                 get<0>(angle.first)->idatm_type(),
                                                 __ring_size(get<0>(angle.first)),
                                                 bin
                                               );
                        }
                        
                        if ( ! bac.count(ab) )
                            bac[ab] = 1;
                        else
                            ++bac[ab];
                 }
                 ba.clear();
        }

        void binDihedrals( const double dihedral_bin_size ) {
                 for ( const auto& dihedral : bd ) {
                        const int bin = floor(dihedral.second / dihedral_bin_size);

                        DihedralBin f=make_tuple(get<0>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<0>(dihedral.first)),
                                                 substitutions.at(get<0>(dihedral.first)),
                                                 get<1>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<1>(dihedral.first)),
                                                 substitutions.at(get<1>(dihedral.first)),
                                                 get<2>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<2>(dihedral.first)),
                                                 substitutions.at(get<2>(dihedral.first)),
                                                 get<3>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<3>(dihedral.first)),
                                                 substitutions.at(get<3>(dihedral.first)),
                                                 bin
                                               );

                        DihedralBin b=make_tuple(get<3>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<3>(dihedral.first)),
                                                 substitutions.at(get<3>(dihedral.first)),
                                                 get<2>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<2>(dihedral.first)),
                                                 substitutions.at(get<2>(dihedral.first)),
                                                 get<1>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<1>(dihedral.first)),
                                                 substitutions.at(get<1>(dihedral.first)),
                                                 get<0>(dihedral.first)->idatm_type(),
                                                 __ring_size(get<0>(dihedral.first)),
                                                 substitutions.at(get<0>(dihedral.first)),
                                                 bin
                                               );
                        DihedralBin db;
                        
                        if ( __comp_atoms(get<0>(dihedral.first), get<3>(dihedral.first) ) < 0 ) {
                                db = f;
                        } else if ( __comp_atoms(get<0>(dihedral.first), get<3>(dihedral.first) ) > 0 ) {
                                db = b;
                        } else if (__comp_atoms(get<1>(dihedral.first), get<2>(dihedral.first) ) < 0) {
                                db = f;
                        } else {
                                db = b;
                        }
                        
                        if ( ! bdc.count(db) )
                            bdc.emplace(db, 1);
                        else
                            ++bdc[db];
                 }
                 bd.clear();
        }

        void binImpropers( const double improper_bin_size ) {
                 for ( const auto& improper : bi ) {
                        const int bin = floor(improper.second / improper_bin_size);

                        const Molib::Atom *atom1, *atom2, *atom3, *atom4;
                        std::tie(atom1, atom2, atom3, atom4) = improper.first;

                        auto thing = [this](const Molib::Atom* atm1, const Molib::Atom* atm2) {
                                return this->__comp_atoms(atm1, atm2);
                        };

                        const Molib::Atom* highest = max( max( atom1, atom3, thing), atom4, thing);
                        const Molib::Atom* lowest  = min( min( atom1, atom3, thing), atom4, thing);
                        const Molib::Atom* middle  = max( min( atom1, atom3, thing),
                                                          min( max(atom1,atom3, thing), atom4, thing),
                                                          thing);

                        DihedralBin ib=make_tuple(lowest->idatm_type(),
                                                 __ring_size(lowest),
                                                 substitutions.at(lowest),
                                                 middle->idatm_type(),
                                                 __ring_size(middle),
                                                 substitutions.at(middle),
                                                 atom2->idatm_type(),
                                                 __ring_size(atom2),
                                                 substitutions.at(atom2),
                                                 highest->idatm_type(),
                                                 __ring_size(highest),
                                                 substitutions.at(highest),
                                                 bin
                                               );

                        if ( ! bic.count(ib) )
                            bic[ib] = 1;
                        else
                            ++bic[ib];
                 }
                 bi.clear();
        }

        void printStretches (std::ostream& os) const {
                for ( const auto& bond : bs ) {
                        __print_atom(get<0>(bond.first), os, ' ');
                        __print_atom(get<1>(bond.first), os, ' ');
                        os << bond.second;
                        os << "\n";
                }
        }

        void printAngles (std::ostream& os) const {
                for ( const auto& angle : ba ) {
                        __print_atom(get<0>(angle.first), os, ' ');
                        __print_atom(get<1>(angle.first), os, ' ');
                        __print_atom(get<2>(angle.first), os, ' ');
                        os << Geom3D::degrees(angle.second);
                        os << "\n";
                }
        }
        
        void printDihedrals (std::ostream& os) const {
                for ( const auto& dihedral : bd ) {
                        __print_atom(get<0>(dihedral.first), os, ' ');
                        __print_atom(get<1>(dihedral.first), os, ' ');
                        __print_atom(get<2>(dihedral.first), os, ' ');
                        __print_atom(get<3>(dihedral.first), os, ' ');
                        os << Geom3D::degrees(dihedral.second);
                        os << "\n";
                }
        }

        void printImpropers (std::ostream& os) const {
                for ( const auto& improper : bi ) {
                        __print_atom(get<0>(improper.first), os, ' ');
                        __print_atom(get<1>(improper.first), os, ' ');
                        __print_atom(get<2>(improper.first), os, ' ');
                        __print_atom(get<3>(improper.first), os, ' ');
                        os << Geom3D::degrees(improper.second);
                        os << "\n";
                }
        }

        void printStretchBins (std::ostream& os) const {
                for ( const auto& bond : bsc ) {
                        os << help::idatm_unmask[get<0>(bond.first)] << " ";
                        os << help::idatm_unmask[get<1>(bond.first)] << " ";
                        os << get<2>(bond.first) << " ";
                        os << bond.second;
                        os << "\n";
                }

        }

        void printAngleBins (std::ostream& os) const {
                for ( const auto& angle : bac ) {
                        os << help::idatm_unmask[get<0>(angle.first)] << " ";
                        os << get<1>(angle.first) << " ";
                        os << help::idatm_unmask[get<2>(angle.first)] << " ";
                        os << get<3>(angle.first) << " ";
                        os << help::idatm_unmask[get<4>(angle.first)] << " ";
                        os << get<5>(angle.first) << " ";
                        os << get<6>(angle.first) << " ";
                        os << angle.second;
                        os << "\n";
                }

        }

        void printDihedralBins (std::ostream& os) const {
                for ( const auto& dihedral : bdc ) {
                        os << help::idatm_unmask[get<0>(dihedral.first)] << " ";
                        os << get<1>(dihedral.first) << " ";
                        os << get<2>(dihedral.first) << " ";
                        os << help::idatm_unmask[get<3>(dihedral.first)] << " ";
                        os << get<4>(dihedral.first) << " ";
                        os << get<5>(dihedral.first) << " ";
                        os << help::idatm_unmask[get<6>(dihedral.first)] << " ";
                        os << get<7>(dihedral.first) << " ";
                        os << get<8>(dihedral.first) << " ";
                        os << help::idatm_unmask[get<9>(dihedral.first)] << " ";
                        os << get<10>(dihedral.first) << " ";
                        os << get<11>(dihedral.first) << " ";
                        os << get<12>(dihedral.first) << " ";
                        os << dihedral.second;
                        os << "\n";
                }

        }
        
        void printImproperBins (std::ostream& os) const {
                for ( const auto& dihedral : bic ) {
                        os << help::idatm_unmask[get<0>(dihedral.first)] << " ";
                        os << get<1>(dihedral.first) << " ";
                        os << get<2>(dihedral.first) << " ";
                        os << help::idatm_unmask[get<3>(dihedral.first)] << " ";
                        os << get<4>(dihedral.first) << " ";
                        os << get<5>(dihedral.first) << " ";
                        os << help::idatm_unmask[get<6>(dihedral.first)] << " ";
                        os << get<7>(dihedral.first) << " ";
                        os << get<8>(dihedral.first) << " ";
                        os << help::idatm_unmask[get<9>(dihedral.first)] << " ";
                        os << get<10>(dihedral.first) << " ";
                        os << get<11>(dihedral.first) << " ";
                        os << get<12>(dihedral.first) << " ";
                        os << dihedral.second;
                        os << "\n";
                }

        }


};

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                MoleculeBondExtractor MBE;

                boost::program_options::variables_map vm;
                po::options_description cmdln_options;

                po::options_description generic;
                generic.add_options()
                ("help,h", "Show this help");

                po::options_description extract_file_options ("Molecular Bond Extraction Files");
                extract_file_options.add_options()
                ("stretch_extract_file", po::value<std::string> ()->default_value("stretches.tbl"),
                 "File to save stretches")
                ("angle_extract_file", po::value<std::string> ()->default_value("angles.tbl"),
                 "File to save angles")
                ("dihedral_extract_file", po::value<std::string> ()->default_value("dihedrals.tbl"),
                 "File to save dihedrals")
                ("improper_extract_file", po::value<std::string> ()->default_value("impropers.tbl"),
                 "File to save impropers")
                ;

                po::options_description binning_options ("Bond Bin Files");
                binning_options.add_options()
                ("stretch_bin_size,s", po::value<double> ()->default_value(0.01),
                 "Stretch bin size")
                ("angle_bin_size,a", po::value<double> ()->default_value(0.01),
                 "Angle bin size")
                ("dihedral_bin_size,d", po::value<double> ()->default_value(0.01),
                 "Dihedral bin size")
                ("improper_bin_size,i", po::value<double> ()->default_value(0.01),
                 "Improper bin size")
                ;

                po::options_description bin_file_options ("Bond Bin Sizes");
                bin_file_options.add_options()
                ("stretch_bin_out", po::value<std::string> ()->default_value("stretch_bins.tbl"),
                 "File to save stretches")
                ("angle_bin_out", po::value<std::string> ()->default_value("angle_bins.tbl"),
                 "File to save angles")
                ("dihedral_bin_out", po::value<std::string> ()->default_value("dihedral_bins.tbl"),
                 "File to save dihedrals")
                ("improper_bin_out", po::value<std::string> ()->default_value("improper_bins.tbl"),
                 "File to save impropers")
                ;

                cmdln_options.add(generic);
                cmdln_options.add(extract_file_options);
                cmdln_options.add(binning_options);
                cmdln_options.add(bin_file_options);

                po::store (po::parse_command_line (argc, argv, cmdln_options), vm);
                po::notify (vm);

                if (vm.count ("help")) {
                        std::cout << cmdln_options << std::endl;
                        return 0;
                }

                std::ofstream stretch_file  (vm["stretch_extract_file"].as<std::string>(),  std::ios::app);
                std::ofstream angle_file    (vm["angle_extract_file"].as<std::string>(),    std::ios::app);
                std::ofstream dihedral_file (vm["dihedral_extract_file"].as<std::string>(), std::ios::app);
                std::ofstream improper_file (vm["improper_extract_file"].as<std::string>(), std::ios::app);

                double stretch_bin = vm["stretch_bin_size"].as<double>();
                double angle_bin = vm["angle_bin_size"].as<double>();
                double dihedral_bin = vm["dihedral_bin_size"].as<double>();
                double improper_bin = vm["improper_bin_size"].as<double>();

                for (int i = 1; i < argc; ++i) {
                        Parser::FileParser input(argv[i], Parser::pdb_read_options::all_models);
                        Molib::Molecules mols;
                        input.parse_molecule(mols);
                        
                        mols.compute_idatm_type();
                        mols.erase_hydrogen();

                        for ( const auto& mol : mols ) {
                                if (!MBE.addMolecule(mol))
                                    continue;

                                MBE.printStretches(stretch_file);
                                MBE.printAngles(stretch_file);
                                MBE.printDihedrals(dihedral_file);
                                MBE.printImpropers(improper_file);

                                MBE.binStretches(stretch_bin);
                                MBE.binAngles(angle_bin);
                                MBE.binDihedrals(dihedral_bin);
                                MBE.binImpropers(improper_bin);
                        }
                }

                stretch_file.close();
                angle_file.close();
                dihedral_file.close();
                improper_file.close();

                std::ofstream stretch_bin_file  (vm["stretch_bin_out"].as<std::string>());
                std::ofstream angle_bin_file    (vm["angle_bin_out"].as<std::string>());
                std::ofstream dihedral_bin_file (vm["dihedral_bin_out"].as<std::string>());
                std::ofstream improper_bin_file (vm["improper_bin_out"].as<std::string>());

                MBE.printStretchBins(stretch_bin_file);
                MBE.printAngleBins(angle_bin_file);
                MBE.printDihedralBins(dihedral_bin_file);
                MBE.printImproperBins(improper_bin_file);

        } catch (po::error& e) {
                std::cerr << "error: " << e.what() << std::endl;
                std::cerr << "Please see the help (-h) for more information" << std::endl;
        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}
