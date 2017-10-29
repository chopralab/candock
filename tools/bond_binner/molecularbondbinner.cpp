#include "molecularbondbinner.hpp"

#include <cmath>
#include "molecularbondextractor.hpp"

using namespace std;

namespace AtomInfo {

void MolecularBondBinner::binStretches(const VBondStretches& bs) {
    for ( const auto& bond : bs ) {
        const size_t bin = std::floor(bond.second / __stretch_bin_size);

        StretchBin sb = make_pair(
            bond.first,
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
            angle.first,
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
            dihedral.first,
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
            improper.first,
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

void MolecularBondBinner::addStretchExtracts(std::istream& is) {
    while (!is.eof()) {
        std::string atom1, atom2;
        int atom_ring_1, atom_ring_2; //FIXME: this type may change
        size_t atom_size_1, atom_size_2;

        double stretch;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;

        is >> stretch;

        if (is.eof())
            break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};

        size_t bin = std::floor( stretch / __stretch_bin_size );

        StretchBin sb = make_pair(
            make_tuple(a1, a2),
            bin);

        if ( ! bsc.count(sb) )
            bsc[sb] = 1;
        else
            ++bsc[sb];

    }
}

void MolecularBondBinner::addAngleExtracts(std::istream& is) {
    while (!is.eof()) {
        std::string atom1, atom2, atom3;
        int atom_ring_1, atom_ring_2, atom_ring_3; //FIXME: this type may change
        size_t atom_size_1, atom_size_2, atom_size_3;

        double angle;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;
        is >> atom3 >> atom_ring_3 >> atom_size_3;

        is >> angle;

        if (is.eof())
            break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};
        atom_info a3 = {help::idatm_mask.at(atom3), atom_ring_3, atom_size_3};

        size_t bin = std::floor( angle / __angle_bin_size );

        AngleBin ab = make_pair(
            make_tuple(a1, a2, a3),
            bin);

        if ( ! bac.count(ab) )
            bac[ab] = 1;
        else
            ++bac[ab];

    }
}

void MolecularBondBinner::addDihedralExtracts(std::istream& is) {
    while (!is.eof()) {
        std::string atom1, atom2, atom3, atom4;
        int atom_ring_1, atom_ring_2, atom_ring_3, atom_ring_4; //FIXME: this type may change
        size_t atom_size_1, atom_size_2, atom_size_3, atom_size_4;

        double dihedral;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;
        is >> atom3 >> atom_ring_3 >> atom_size_3;
        is >> atom4 >> atom_ring_4 >> atom_size_4;

        is >> dihedral;

        if (is.eof())
            break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};
        atom_info a3 = {help::idatm_mask.at(atom3), atom_ring_3, atom_size_3};
        atom_info a4 = {help::idatm_mask.at(atom4), atom_ring_4, atom_size_4};

        int bin = std::floor( dihedral / __dihedral_bin_size );

        DihedralBin db = make_pair(
            make_tuple(a1, a2, a3, a4),
            bin);

        if ( ! bdc.count(db) )
            bdc[db] = 1;
        else
            ++bdc[db];

    }
}

void MolecularBondBinner::addImproperExtracts(std::istream& is) {
    while (!is.eof()) {
        std::string atom1, atom2, atom3, atom4;
        int atom_ring_1, atom_ring_2, atom_ring_3, atom_ring_4; //FIXME: this type may change
        size_t atom_size_1, atom_size_2, atom_size_3, atom_size_4;

        double improper;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;
        is >> atom3 >> atom_ring_3 >> atom_size_3;
        is >> atom4 >> atom_ring_4 >> atom_size_4;

        is >> improper;

        if (is.eof())
            break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};
        atom_info a3 = {help::idatm_mask.at(atom3), atom_ring_3, atom_size_3};
        atom_info a4 = {help::idatm_mask.at(atom4), atom_ring_4, atom_size_4};

        int bin = std::floor( improper / __improper_bin_size );

        DihedralBin ib = make_pair(
            make_tuple(a1, a2, a3, a4),
            bin);

        if ( ! bic.count(ib) )
            bic[ib] = 1;
        else
            ++bic[ib];

    }
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
    for ( const auto& improper : bic ) {
        os << get<0>(improper.first.first) << " ";
        os << get<1>(improper.first.first) << " ";
        os << get<2>(improper.first.first) << " ";
        os << get<3>(improper.first.first) << " ";
        os << improper.first.second << " ";
        os << improper.second;
        os << "\n";
    }

}

}
