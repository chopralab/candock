#ifndef MOLECULARBONDEXTRACTOR_HPP
#define MOLECULARBONDEXTRACTOR_HPP

#include <tuple>
#include <map>
#include "molib/atom.hpp"
#include "molib/molecule.hpp"

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
    // FIXME: consider changing the types here

    struct bond_info {
        int idatm_type;
        int ring_size;
        size_t substitutions;
        
        friend bool operator< (const bond_info& lhs, const bond_info& rhs) {
            if ( lhs.idatm_type < rhs.idatm_type ) {
                return -4;
            } else if ( lhs.idatm_type > rhs.idatm_type ) {
                return 4;
            }

            if ( lhs.ring_size < rhs.ring_size ) {
                return -3;
            } else if ( lhs.ring_size > rhs.ring_size) {
                return 3;
            }

            if ( lhs.substitutions < rhs.substitutions ) {
                return -2;
            } else if (lhs.substitutions > rhs.substitutions) {
                return 2;
            }

            return 0;
        }
    };

    typedef std::tuple<bond_info, bond_info, size_t> StretchBin;
    typedef std::map< StretchBin, size_t >           StretchCounts;
    typedef std::tuple<bond_info, bond_info,
                       bond_info, size_t >           AngleBin;
    typedef std::map< AngleBin, size_t >             AngleCounts;
    typedef std::tuple<bond_info, bond_info,
                       bond_info, bond_info, size_t> DihedralBin;
    typedef std::map< DihedralBin, size_t >          DihedralCounts;

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

    int __comp_atoms ( const Molib::Atom* atom1, const Molib::Atom* atom2) const;

    int __ring_size( const Molib::Atom* atom ) const;

    void __print_atom( const Molib::Atom* atom, std::ostream& os, char delim ) const;

    bond_info __make_bond_info( const Molib::Atom* atom);
    

public:

    void reset();

    void addStretch ( BondStretch stretch_atoms );
    void addAngle ( BondAngle angle_atoms);
    void addDihedral (BondDihedral dihedral_atoms);
    void addImproper (BondDihedral improper_atoms);
    bool addMolecule( const Molib::Molecule& mol );

    void binStretches( const double bond_bin_size );
    void binAngles( const double angle_bin_size );
    void binDihedrals( const double dihedral_bin_size );
    void binImpropers( const double improper_bin_size );

    void printStretches (std::ostream& os) const;
    void printAngles (std::ostream& os) const;
    void printDihedrals (std::ostream& os) const;
    void printImpropers (std::ostream& os) const;

    void printStretchBins (std::ostream& os) const;
    void printAngleBins (std::ostream& os) const;
    void printDihedralBins (std::ostream& os) const;
    void printImproperBins (std::ostream& os) const;

};

#endif
