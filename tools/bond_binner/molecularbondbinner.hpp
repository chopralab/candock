#ifndef MOLECULARBONDBINNER_HPP
#define MOLECULARBONDBINNER_HPP

#include <ostream>
#include <vector>
#include "atominfotypes.hpp"

namespace AtomInfo {

class MolecularBondExtractor;

class MolecularBondBinner
{
    double __stretch_bin_size;
    double __angle_bin_size;
    double __dihedral_bin_size;
    double __improper_bin_size;

    void __addStretch( const BondStretchBin& sb, size_t inc = 1);
    void __addAngle( const BondAngleBin& sb, size_t inc = 1);
    void __addDihedral( const BondDihedralBin& sb, size_t inc = 1);
    void __addImproper( const BondDihedralBin& sb, size_t inc = 1);

    void binStretches(const VBondStretchValues& bs);
    void binAngles( const VBondAngleValues& ba );
    void binDihedrals( const VBondDihedralValues& bd );
    void binImpropers( const VBondDihedralValues& bi );

    StretchCounts  bsc;
    AngleCounts    bac;
    DihedralCounts bdc;
    DihedralCounts bic;

public:

    MolecularBondBinner(double stretch_bin_size, double angle_bin_size,
                        double dihedral_bin_size, double improper_bin_size) : 
                        __stretch_bin_size(stretch_bin_size),
                        __angle_bin_size(angle_bin_size),
                        __dihedral_bin_size(dihedral_bin_size),
                        __improper_bin_size(improper_bin_size){}

    void addExtracts( const MolecularBondExtractor& mbe);
    void addExtracts( const MolecularBondBinner& mbb);

    void addStretchExtracts( std::istream& os, bool is_bin_file = false);
    void addAngleExtracts( std::istream& os,   bool is_bin_file = false);
    void addDihedralExtracts( std::istream& os,bool is_bin_file = false);
    void addImproperExtracts( std::istream& os,bool is_bin_file = false);

    void printStretchBins (std::ostream& os) const;
    void printAngleBins (std::ostream& os) const;
    void printDihedralBins (std::ostream& os) const;
    void printImproperBins (std::ostream& os) const;

};

}

#endif // MOLECULARBONDBINNER_HPP
