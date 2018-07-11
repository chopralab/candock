#include "molecularbondbinner.hpp"

#include <cmath>
#include "molecularbondextractor.hpp"

using namespace std;
using namespace candock;

namespace AtomInfo {

void MolecularBondBinner::__addStretch(const BondStretchBin& sb, size_t inc) {
    if (!bsc.count(sb))
        bsc[sb] = inc;
    else
        bsc[sb] += inc;
}

void MolecularBondBinner::__addAngle(const BondAngleBin& sb, size_t inc) {
    if (!bac.count(sb))
        bac[sb] = inc;
    else
        bac[sb] += inc;
}

void MolecularBondBinner::__addDihedral(const BondDihedralBin& sb, size_t inc) {
    if (!bdc.count(sb))
        bdc[sb] = inc;
    else
        bdc[sb] += inc;
}

void MolecularBondBinner::__addImproper(const BondDihedralBin& sb, size_t inc) {
    if (!bic.count(sb))
        bic[sb] = inc;
    else
        bic[sb] += inc;
}

void MolecularBondBinner::binStretches(const VBondStretchValues& bs) {
    for (const auto& bond : bs) {
        const size_t bin = std::lround(bond.second / __stretch_bin_size);

        BondStretchBin sb = make_pair(bond.first, bin);

        __addStretch(sb);
    }
}

void MolecularBondBinner::binAngles(const VBondAngleValues& ba) {
    for (const auto& angle : ba) {
        double angle_val = angle.second;

        while (angle_val < 0) angle_val += M_PI * 2;

        const size_t bin = std::lround(angle_val / __angle_bin_size);

        BondAngleBin ab = make_pair(angle.first, bin);

        __addAngle(ab);
    }
}

void MolecularBondBinner::binDihedrals(const VBondDihedralValues& bd) {
    for (const auto& dihedral : bd) {
        // Dihedrals can be negative, fix it!
        double dihedral_val = dihedral.second;

        while (dihedral_val < 0) dihedral_val += M_PI * 2;

        const size_t bin = std::lround(dihedral_val / __dihedral_bin_size);

        BondDihedralBin db = make_pair(dihedral.first, bin);

        __addDihedral(db);
    }
}

void MolecularBondBinner::binImpropers(const VBondDihedralValues& bi) {
    for (const auto& improper : bi) {
        double improper_val = improper.second;

        while (improper_val < 0) improper_val += M_PI * 2;

        const size_t bin = std::lround(improper_val / __improper_bin_size);

        BondDihedralBin ib = make_pair(improper.first, bin);

        __addImproper(ib);
    }
}

void MolecularBondBinner::addExtracts(const MolecularBondExtractor& mbe) {
    binStretches(mbe.stretches());
    binAngles(mbe.angles());
    binDihedrals(mbe.dihedrals());
    binImpropers(mbe.impropers());
}

void MolecularBondBinner::addExtracts(const MolecularBondBinner& mbb) {
    for (const auto& kv : mbb.bsc) {
        __addStretch(kv.first, kv.second);
    }

    for (const auto& kv : mbb.bac) {
        __addAngle(kv.first, kv.second);
    }

    for (const auto& kv : mbb.bdc) {
        __addDihedral(kv.first, kv.second);
    }

    for (const auto& kv : mbb.bic) {
        __addImproper(kv.first, kv.second);
    }
}

void MolecularBondBinner::addStretchExtracts(std::istream& is,
                                             bool is_bin_file) {
    while (!is.eof()) {
        std::string atom1, atom2;
        size_t atom_ring_1, atom_ring_2;
        size_t atom_size_1, atom_size_2;

        double stretch;
        size_t bin, count = 1;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;

        if (is_bin_file)
            is >> bin >> count;
        else
            is >> stretch;

        if (is.eof()) break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};

        if (!is_bin_file) bin = std::lround(stretch / __stretch_bin_size);

        BondStretchBin sb = make_pair(make_tuple(a1, a2), bin);

        __addStretch(sb, count);
    }
}

void MolecularBondBinner::addAngleExtracts(std::istream& is, bool is_bin_file) {
    while (!is.eof()) {
        std::string atom1, atom2, atom3;
        size_t atom_ring_1, atom_ring_2, atom_ring_3;
        size_t atom_size_1, atom_size_2, atom_size_3;

        double angle;
        size_t bin, count = 1;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;
        is >> atom3 >> atom_ring_3 >> atom_size_3;

        if (is_bin_file)
            is >> bin >> count;
        else
            is >> angle;

        if (is.eof()) break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};
        atom_info a3 = {help::idatm_mask.at(atom3), atom_ring_3, atom_size_3};

        if (!is_bin_file) {
            // Just in case...
            while (angle < 0) angle += M_PI * 2;

            bin = std::lround(angle / __angle_bin_size);
        }

        BondAngleBin ab = make_pair(make_tuple(a1, a2, a3), bin);

        __addAngle(ab, count);
    }
}

void MolecularBondBinner::addDihedralExtracts(std::istream& is,
                                              bool is_bin_file) {
    while (!is.eof()) {
        std::string atom1, atom2, atom3, atom4;
        size_t atom_ring_1, atom_ring_2, atom_ring_3, atom_ring_4;
        size_t atom_size_1, atom_size_2, atom_size_3, atom_size_4;

        double dihedral;
        size_t bin, count = 1;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;
        is >> atom3 >> atom_ring_3 >> atom_size_3;
        is >> atom4 >> atom_ring_4 >> atom_size_4;

        if (is_bin_file)
            is >> bin >> count;
        else
            is >> dihedral;

        if (is.eof()) break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};
        atom_info a3 = {help::idatm_mask.at(atom3), atom_ring_3, atom_size_3};
        atom_info a4 = {help::idatm_mask.at(atom4), atom_ring_4, atom_size_4};

        if (!is_bin_file) {
            // Adjust for negative dihedrals
            while (dihedral < 0) dihedral += M_PI * 2;

            bin = std::lround(dihedral / __dihedral_bin_size);
        }

        BondDihedralBin db = make_pair(make_tuple(a1, a2, a3, a4), bin);

        __addDihedral(db, count);
    }
}

void MolecularBondBinner::addImproperExtracts(std::istream& is,
                                              bool is_bin_file) {
    while (!is.eof()) {
        std::string atom1, atom2, atom3, atom4;
        size_t atom_ring_1, atom_ring_2, atom_ring_3, atom_ring_4;
        size_t atom_size_1, atom_size_2, atom_size_3, atom_size_4;

        double improper;
        size_t bin, count = 1;

        is >> atom1 >> atom_ring_1 >> atom_size_1;
        is >> atom2 >> atom_ring_2 >> atom_size_2;
        is >> atom3 >> atom_ring_3 >> atom_size_3;
        is >> atom4 >> atom_ring_4 >> atom_size_4;

        if (is_bin_file)
            is >> bin >> count;
        else
            is >> improper;
        ;

        if (is.eof()) break;

        atom_info a1 = {help::idatm_mask.at(atom1), atom_ring_1, atom_size_1};
        atom_info a2 = {help::idatm_mask.at(atom2), atom_ring_2, atom_size_2};
        atom_info a3 = {help::idatm_mask.at(atom3), atom_ring_3, atom_size_3};
        atom_info a4 = {help::idatm_mask.at(atom4), atom_ring_4, atom_size_4};

        if (!is_bin_file) {
            while (improper < 0) improper += M_PI * 2;

            bin = std::lround(improper / __improper_bin_size);
        }

        BondDihedralBin ib = make_pair(make_tuple(a1, a2, a3, a4), bin);

        __addImproper(ib, count);
    }
}

void MolecularBondBinner::printStretchBins(std::ostream& os) const {
    for (const auto& bond : bsc) {
        os << get<0>(bond.first.first) << " ";
        os << get<1>(bond.first.first) << " ";
        os << bond.first.second << " ";
        os << bond.second;
        os << "\n";
    }
}

void MolecularBondBinner::printAngleBins(std::ostream& os) const {
    for (const auto& angle : bac) {
        os << get<0>(angle.first.first) << " ";
        os << get<1>(angle.first.first) << " ";
        os << get<2>(angle.first.first) << " ";
        os << angle.first.second << " ";
        os << angle.second;
        os << "\n";
    }
}

void MolecularBondBinner::printDihedralBins(std::ostream& os) const {
    for (const auto& dihedral : bdc) {
        os << get<0>(dihedral.first.first) << " ";
        os << get<1>(dihedral.first.first) << " ";
        os << get<2>(dihedral.first.first) << " ";
        os << get<3>(dihedral.first.first) << " ";
        os << dihedral.first.second << " ";
        os << dihedral.second;
        os << "\n";
    }
}

void MolecularBondBinner::printImproperBins(std::ostream& os) const {
    for (const auto& improper : bic) {
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
