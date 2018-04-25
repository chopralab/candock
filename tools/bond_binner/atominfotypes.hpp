#ifndef ATOMINFOTYPES_HPP
#define ATOMINFOTYPES_HPP

#include <cstddef>
#include <ostream>
#include <map>
#include <tuple>
#include "helper/help.hpp"

using std::size_t;

namespace AtomInfo {

    // typedefs for binned Bond information
    // NOTE: Bondorder not included (yet?)
    // FIXME: consider changing the types here

    struct atom_info {
        int idatm_type;
        size_t ring_size;
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

        friend std::ostream& operator<< (std::ostream& os, const atom_info& ai) {
            os << std::setw(3) << help::idatm_unmask[ ai.idatm_type ] << ' ';
            os << ai.ring_size << ' ';
            os << ai.substitutions << ' ';
            return os;
        }
    };

    // typedefs for raw Bond information
    typedef std::tuple< atom_info, atom_info>    BondStretchAtoms;
    typedef std::pair< BondStretchAtoms, double> BondStretchValue;
    typedef std::vector<BondStretchValue>        VBondStretchValues;

    typedef std::tuple< atom_info, atom_info,
                                   atom_info>  BondAngleAtoms;
    typedef std::pair< BondAngleAtoms, double> BondAngleValue;
    typedef std::vector<BondAngleValue>        VBondAngleValues;
    
    typedef std::tuple< atom_info, atom_info,
                        atom_info, atom_info>     BondDihedralAtoms;
    typedef std::pair< BondDihedralAtoms, double> BondDihedralValue;
    typedef std::vector<BondDihedralValue>        VBondDihedralValues;

    // typedefs for binned data
    typedef std::pair<BondStretchAtoms, size_t>  BondStretchBin;
    typedef std::map< BondStretchBin, size_t >   StretchCounts;

    typedef std::pair<BondAngleAtoms, size_t >   BondAngleBin;
    typedef std::map< BondAngleBin, size_t >     AngleCounts;

    typedef std::pair<BondDihedralAtoms, size_t> BondDihedralBin;
    typedef std::map< BondDihedralBin, size_t >  DihedralCounts;
}

#endif
