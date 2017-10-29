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

    void binStretches(const VBondStretches& bs);
    void binAngles( const VBondAngles& ba );
    void binDihedrals( const VBondDihedrals& bd );
    void binImpropers( const VBondDihedrals& bi );

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
    void addStretchExtracts( std::istream& os);
    void addAngleExtracts( std::istream& os);
    void addDihedralExtracts( std::istream& os);
    void addImproperExtracts( std::istream& os);

    void printStretchBins (std::ostream& os) const;
    void printAngleBins (std::ostream& os) const;
    void printDihedralBins (std::ostream& os) const;
    void printImproperBins (std::ostream& os) const;

};

}

#endif // MOLECULARBONDBINNER_HPP
