#include "molecularbondbinner.hpp"

#include <cmath>
#include "molecularbondextractor.hpp"

using namespace std;

namespace AtomInfo {

void MolecularBondBinner::binStretches(const VBondStretches& bs) {
    for ( const auto& bond : bs ) {
        const size_t bin = std::floor(bond.second / __stretch_bin_size);

        StretchBin sb = make_pair(
            make_tuple(
                get<0>(bond.first),
                get<1>(bond.first)
            ),
            bin);

        if ( ! bsc.count(sb) )
            bsc[sb] = 1;
        else
            ++bsc[sb];
    }
}

void MolecularBondBinner::binAngles(const VBondAngles& ba) {
    for ( const auto& angle : ba ) {
        const size_t bin = std::floor(angle.second / __angle_bin_size);

        AngleBin ab = make_pair(
            make_tuple(
                get<0>(angle.first),
                get<1>(angle.first),
                get<2>(angle.first)
            ),
            bin);

        if ( ! bac.count(ab) )
            bac[ab] = 1;
        else
            ++bac[ab];
    }
}

void MolecularBondBinner::binDihedrals(const VBondDihedrals& bd) {

    for ( const auto& dihedral : bd ) {
        const int bin = std::floor(dihedral.second / __dihedral_bin_size);

        DihedralBin db = make_pair(
            make_tuple(
                get<0>(dihedral.first),
                get<1>(dihedral.first),
                get<2>(dihedral.first),
                get<3>(dihedral.first)
            ),
            bin);

        if ( ! bdc.count(db) )
            bdc.emplace(db, 1);
        else
            ++bdc[db];
    }
}

void MolecularBondBinner::binImpropers(const VBondDihedrals& bi) {

    for ( const auto& improper : bi ) {
        const int bin = std::floor(improper.second / __improper_bin_size);

        DihedralBin ib = make_pair(
            make_tuple(
                get<0>(improper.first),
                get<1>(improper.first),
                get<2>(improper.first),
                get<3>(improper.first)
            ),
            bin);

        if ( ! bic.count(ib) )
            bic[ib] = 1;
        else
            ++bic[ib];
    }
}

void MolecularBondBinner::addExtracts( const MolecularBondExtractor& mbe) {
    binStretches(mbe.stretches());
    binAngles(mbe.angles());
    binDihedrals(mbe.dihedrals());
    binImpropers(mbe.impropers());
}


void MolecularBondBinner::printStretchBins (std::ostream& os) const {
    for ( const auto& bond : bsc ) {
        os << get<0>(bond.first.first) << " ";
        os << get<1>(bond.first.first) << " ";
        os << bond.first.second << " ";
        os << bond.second;
        os << "\n";
    }

}

void MolecularBondBinner::printAngleBins (std::ostream& os) const {
    for ( const auto& angle : bac ) {
        os << get<0>(angle.first.first) << " ";
        os << get<1>(angle.first.first) << " ";
        os << get<2>(angle.first.first) << " ";
        os << angle.first.second << " ";
        os << angle.second;
        os << "\n";
    }

}

void MolecularBondBinner::printDihedralBins (std::ostream& os) const {
    for ( const auto& dihedral : bdc ) {
        os << get<0>(dihedral.first.first) << " ";
        os << get<1>(dihedral.first.first) << " ";
        os << get<2>(dihedral.first.first) << " ";
        os << get<3>(dihedral.first.first) << " ";
        os << dihedral.first.second << " ";
        os << dihedral.second;
        os << "\n";
    }

}

void MolecularBondBinner::printImproperBins (std::ostream& os) const {
    for ( const auto& dihedral : bic ) {
        os << get<0>(dihedral.first.first) << " ";
        os << get<1>(dihedral.first.first) << " ";
        os << get<2>(dihedral.first.first) << " ";
        os << get<3>(dihedral.first.first) << " ";
        os << dihedral.first.second << " ";
        os << dihedral.second;
        os << "\n";
    }

}

}
