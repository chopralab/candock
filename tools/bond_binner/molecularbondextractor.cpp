#include "molecularbondextractor.hpp"


int MoleculeBondExtractor::__comp_atoms ( const Molib::Atom* atom1,
                                          const Molib::Atom* atom2) const {

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

int MoleculeBondExtractor::__ring_size( const Molib::Atom* atom ) const {
    const std::vector<size_t>& ring_vec = all_rings_sizes.at(atom);

    // Not in a ring
    if (ring_vec.size() == 0) {
        return 0;
    }

    // Bridged or Spiro
    if (ring_vec.size() > 1) {

        // Spiro must be carbon
        if ( atom->element() != Molib::Element::C)
            return 2;

        // Spiro must have four non-hydrogen neighbors
        if ( substitutions.at(atom) != 4 ) {
            return 2;
        }

        // "Difficult" case, we need to check the neighbors
        for ( const Molib::Bond* bond : atom->get_bonds() ) {
            const Molib::Atom& other = bond->second_atom(*atom);

            // Spiro must be connected to only rings
            if ( all_rings_sizes.at(&other).size() == 0 ) {
                return 2;
            }

            // Vast majority of cases present here,
            // neighbor is only in one ring
            if ( all_rings_sizes.at(&other).size() == 1 ) {
                continue;
            }

            // We need to see how many rings the two share
            size_t shared_rings = 0;
            for ( const auto& thing : all_rings ) {

                // Are they in the same ring?
                if ( thing.count( const_cast<Molib::Atom*>(atom) ) &&
                        thing.count( const_cast<Molib::Atom*>(&other) ) ) {
                    shared_rings++;
                }
            }

            // FIXME: should be adjusted for large ring structures (by subtracting total rings?).
            if ( shared_rings >= 2 )
                return 2;
        }

        return 1;
    }

    return ring_vec.at(0);
}

void MoleculeBondExtractor::__print_atom( const Molib::Atom* atom, std::ostream& os, char delim ) const {
    os << help::idatm_unmask[ atom->idatm_type() ] << delim;
    os << __ring_size(atom) << delim;
    os << substitutions.at(atom) << delim;

}

MoleculeBondExtractor::bond_info MoleculeBondExtractor::__make_bond_info(const Molib::Atom* atom) {
    bond_info a = { atom->idatm_type(),
                    __ring_size(atom),
                    substitutions.at(atom)
    };

    return a;
}

void MoleculeBondExtractor::reset() {
    substitutions.clear();
    all_rings.clear();
    all_rings_sizes.clear();
}

void MoleculeBondExtractor::addStretch ( BondStretch stretch_atoms ) {

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

void MoleculeBondExtractor::addAngle ( BondAngle angle_atoms) {

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

void MoleculeBondExtractor::addDihedral (BondDihedral dihedral_atoms) {

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

void MoleculeBondExtractor::addImproper (BondDihedral improper_atoms) {

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


bool MoleculeBondExtractor::addMolecule( const Molib::Molecule& mol ) {

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

void MoleculeBondExtractor::binStretches( const double bond_bin_size ) {
    for ( const auto& bond : bs ) {
        const size_t bin = floor(bond.second / bond_bin_size);

        const bond_info a0 = __make_bond_info(get<0>(bond.first));
        const bond_info a1 = __make_bond_info(get<1>(bond.first));

        StretchBin sb;
        if ( __comp_atoms(get<0>(bond.first), get<1>(bond.first) ) < 0 ) {
            sb = make_tuple( a0, a1, bin );
        } else {
            sb = make_tuple( a1, a0, bin );
        }

        if ( ! bsc.count(sb) )
            bsc[sb] = 1;
        else
            ++bsc[sb];
    }
    bs.clear();
}

void MoleculeBondExtractor::binAngles( const double angle_bin_size ) {
    for ( const auto& angle : ba ) {
        const size_t bin = floor(angle.second / angle_bin_size);

        const bond_info a0 = __make_bond_info(get<0>(angle.first));
        const bond_info a1 = __make_bond_info(get<1>(angle.first));
        const bond_info a2 = __make_bond_info(get<2>(angle.first));

        AngleBin ab;
        if ( __comp_atoms(get<0>(angle.first), get<2>(angle.first) ) < 0 ) {
            ab = make_tuple( a0, a1, a2, bin );
        } else {
            ab = make_tuple( a2, a1, a0, bin );
        }

        if ( ! bac.count(ab) )
            bac[ab] = 1;
        else
            ++bac[ab];
    }
    ba.clear();
}

void MoleculeBondExtractor::binDihedrals( const double dihedral_bin_size ) {
    for ( const auto& dihedral : bd ) {
        const int bin = floor(dihedral.second / dihedral_bin_size);

        const bond_info a0 = __make_bond_info(get<0>(dihedral.first));
        const bond_info a1 = __make_bond_info(get<1>(dihedral.first));
        const bond_info a2 = __make_bond_info(get<2>(dihedral.first));
        const bond_info a3 = __make_bond_info(get<3>(dihedral.first));


        DihedralBin f = make_tuple(a0,a1,a2,a3,bin);
        DihedralBin b = make_tuple(a3,a2,a1,a0,bin);
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

void MoleculeBondExtractor::binImpropers( const double improper_bin_size ) {
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

        const bond_info a0 = __make_bond_info(lowest);
        const bond_info a1 = __make_bond_info(middle);
        const bond_info a2 = __make_bond_info(atom2);
        const bond_info a3 = __make_bond_info(highest);

        DihedralBin ib=make_tuple(a0, a1, a2, a3, bin);

        if ( ! bic.count(ib) )
            bic[ib] = 1;
        else
            ++bic[ib];
    }
    bi.clear();
}

void MoleculeBondExtractor::printStretches (std::ostream& os) const {
    for ( const auto& bond : bs ) {
        __print_atom(get<0>(bond.first), os, ' ');
        __print_atom(get<1>(bond.first), os, ' ');
        os << bond.second;
        os << "\n";
    }
}

void MoleculeBondExtractor::printAngles (std::ostream& os) const {
    for ( const auto& angle : ba ) {
        __print_atom(get<0>(angle.first), os, ' ');
        __print_atom(get<1>(angle.first), os, ' ');
        __print_atom(get<2>(angle.first), os, ' ');
        os << Geom3D::degrees(angle.second);
        os << "\n";
    }
}

void MoleculeBondExtractor::printDihedrals (std::ostream& os) const {
    for ( const auto& dihedral : bd ) {
        __print_atom(get<0>(dihedral.first), os, ' ');
        __print_atom(get<1>(dihedral.first), os, ' ');
        __print_atom(get<2>(dihedral.first), os, ' ');
        __print_atom(get<3>(dihedral.first), os, ' ');
        os << Geom3D::degrees(dihedral.second);
        os << "\n";
    }
}

void MoleculeBondExtractor::printImpropers (std::ostream& os) const {
    for ( const auto& improper : bi ) {
        __print_atom(get<0>(improper.first), os, ' ');
        __print_atom(get<1>(improper.first), os, ' ');
        __print_atom(get<2>(improper.first), os, ' ');
        __print_atom(get<3>(improper.first), os, ' ');
        os << Geom3D::degrees(improper.second);
        os << "\n";
    }
}

void MoleculeBondExtractor::printStretchBins (std::ostream& os) const {
    for ( const auto& bond : bsc ) {
        os << help::idatm_unmask[get<0>(bond.first).idatm_type] << " ";
        os << help::idatm_unmask[get<1>(bond.first).idatm_type] << " ";
        os << get<2>(bond.first) << " ";
        os << bond.second;
        os << "\n";
    }

}

void MoleculeBondExtractor::printAngleBins (std::ostream& os) const {
    for ( const auto& angle : bac ) {
        os << help::idatm_unmask[get<0>(angle.first).idatm_type] << " ";
        os << get<0>(angle.first).ring_size << " ";
        os << help::idatm_unmask[get<1>(angle.first).idatm_type] << " ";
        os << get<1>(angle.first).ring_size << " ";
        os << help::idatm_unmask[get<2>(angle.first).idatm_type] << " ";
        os << get<2>(angle.first).ring_size << " ";
        os << get<3>(angle.first) << " ";
        os << angle.second;
        os << "\n";
    }

}

void MoleculeBondExtractor::printDihedralBins (std::ostream& os) const {
    for ( const auto& dihedral : bdc ) {
        os << help::idatm_unmask[get<0>(dihedral.first).idatm_type] << " ";
        os << get<0>(dihedral.first).ring_size << " ";
        os << get<0>(dihedral.first).substitutions << " ";
        os << help::idatm_unmask[get<1>(dihedral.first).idatm_type] << " ";
        os << get<1>(dihedral.first).ring_size << " ";
        os << get<1>(dihedral.first).substitutions << " ";
        os << help::idatm_unmask[get<2>(dihedral.first).idatm_type] << " ";
        os << get<2>(dihedral.first).ring_size << " ";
        os << get<2>(dihedral.first).substitutions << " ";
        os << help::idatm_unmask[get<3>(dihedral.first).idatm_type] << " ";
        os << get<3>(dihedral.first).ring_size << " ";
        os << get<3>(dihedral.first).substitutions << " ";
        os << get<4>(dihedral.first) << " ";
        os << dihedral.second;
        os << "\n";
    }

}

void MoleculeBondExtractor::printImproperBins (std::ostream& os) const {
    for ( const auto& dihedral : bic ) {
        os << help::idatm_unmask[get<0>(dihedral.first).idatm_type] << " ";
        os << get<0>(dihedral.first).ring_size << " ";
        os << get<0>(dihedral.first).substitutions << " ";
        os << help::idatm_unmask[get<1>(dihedral.first).idatm_type] << " ";
        os << get<1>(dihedral.first).ring_size << " ";
        os << get<1>(dihedral.first).substitutions << " ";
        os << help::idatm_unmask[get<2>(dihedral.first).idatm_type] << " ";
        os << get<2>(dihedral.first).ring_size << " ";
        os << get<2>(dihedral.first).substitutions << " ";
        os << help::idatm_unmask[get<3>(dihedral.first).idatm_type] << " ";
        os << get<3>(dihedral.first).ring_size << " ";
        os << get<3>(dihedral.first).substitutions << " ";
        os << get<4>(dihedral.first) << " ";
        os << dihedral.second;
        os << "\n";
    }

}
