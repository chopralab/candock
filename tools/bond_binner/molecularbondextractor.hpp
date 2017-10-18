#ifndef MOLECULARBONDEXTRACTOR_HPP
#define MOLECULARBONDEXTRACTOR_HPP

#include <tuple>
#include <map>
#include "molib/atom.hpp"
#include "molib/molecule.hpp"

class MoleculeBondExtractor {

public:

    // typedefs for binned Bond information
    // NOTE: Bondorder not included (yet?)
    // FIXME: consider changing the types here

    struct atom_info {
        int idatm_type;
        int ring_size;
        size_t substitutions;
        
        friend bool operator< (const atom_info& lhs, const atom_info& rhs) {
            if ( lhs.idatm_type < rhs.idatm_type ) {
                return true;
            } else if ( lhs.idatm_type > rhs.idatm_type ) {
                return false;
            }

            if ( lhs.ring_size < rhs.ring_size ) {
                return true;
            } else if ( lhs.ring_size > rhs.ring_size) {
                return false;
            }

            if ( lhs.substitutions < rhs.substitutions ) {
                return true;
            } else if (lhs.substitutions > rhs.substitutions) {
                return false;
            }

            return false;
        }

        friend std::ostream& operator<< (std::ostream& os, const atom_info& ai);
    };

    // typedefs for raw Bond information
    typedef std::tuple< atom_info, atom_info> BondStretch;

    typedef std::pair< BondStretch, double>    BondStretches;

    typedef std::tuple< atom_info, atom_info,
                                   atom_info> BondAngle;
    typedef std::pair< BondAngle, double>     BondAngles;
    
    typedef std::tuple< atom_info, atom_info,
                        atom_info, atom_info> BondDihedral;
    typedef std::pair< BondDihedral, double>  BondDihedrals;

    typedef std::pair<BondStretch, size_t>  StretchBin;
    typedef std::map< StretchBin, size_t >  StretchCounts;
    typedef std::pair<BondAngle, size_t >   AngleBin;
    typedef std::map< AngleBin, size_t >    AngleCounts;
    typedef std::pair<BondDihedral, size_t> DihedralBin;
    typedef std::map< DihedralBin, size_t > DihedralCounts;

private:
    std::vector<BondStretches> bs;
    std::vector<BondAngles>    ba;
    std::vector<BondDihedrals> bd;
    std::vector<BondDihedrals> bi; // Impropers

    StretchCounts  bsc;
    AngleCounts    bac;
    DihedralCounts bdc;
    DihedralCounts bic;

    Glib::Graph<Molib::Atom>::Cycles all_rings;
    Glib::Graph<Molib::Atom>::VertexRingMap all_rings_sizes;
    std::map<const Molib::Atom*, size_t> substitutions;

    int __ring_size( const Molib::Atom* atom ) const;

    void __print_atom( const Molib::Atom* atom, std::ostream& os, char delim ) const;

    atom_info __make_atom_info( const Molib::Atom* atom);
    
    void addStretch (const Molib::Atom* atom1, const Molib::Atom* atom2 );
    void addAngle   (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3);
    void addDihedral(const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3, const Molib::Atom* atom4);
    void addImproper(const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3, const Molib::Atom* atom4);

public:

    void reset();

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
