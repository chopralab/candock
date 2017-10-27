#ifndef ATOMINFOTYPES_HPP
#define ATOMINFOTYPES_HPP

#include <ostream>
#include <map>
#include <tuple>

namespace AtomInfo {

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

    typedef std::pair< BondStretch, double>   BondStretches;
    typedef std::vector<BondStretches>        VBondStretches;

    typedef std::tuple< atom_info, atom_info,
                                   atom_info> BondAngle;
    typedef std::pair< BondAngle, double>     BondAngles;
    typedef std::vector<BondAngles>           VBondAngles;
    
    typedef std::tuple< atom_info, atom_info,
                        atom_info, atom_info> BondDihedral;
    typedef std::pair< BondDihedral, double>  BondDihedrals;
    typedef std::vector<BondDihedrals>        VBondDihedrals;

    typedef std::pair<BondStretch, size_t>  StretchBin;
    typedef std::map< StretchBin, size_t >  StretchCounts;

    typedef std::pair<BondAngle, size_t >   AngleBin;
    typedef std::map< AngleBin, size_t >    AngleCounts;

    // Dihedrals can be negative
    typedef std::pair<BondDihedral, int> DihedralBin;
    typedef std::map< DihedralBin, size_t > DihedralCounts;
}

#endif
