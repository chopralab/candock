#include "molecularbondextractor.hpp"

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

MoleculeBondExtractor::atom_info MoleculeBondExtractor::__make_atom_info(const Molib::Atom* atom) {
    atom_info a = { atom->idatm_type(),
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

void MoleculeBondExtractor::addStretch (const Molib::Atom* atom1, const Molib::Atom* atom2) {

    if (atom1 == atom2)
        return;

    atom_info a0 = __make_atom_info(atom1);
    atom_info a1 = __make_atom_info(atom2);
    
    BondStretch s_new;

    if ( a0 < a1 ) {
        s_new = make_tuple(a0, a1);
    } else {
        s_new = make_tuple(a1, a0);
    }

    double distance = atom1->crd().distance(atom2->crd());
    bs.push_back(make_pair(s_new,distance));
}

void MoleculeBondExtractor::addAngle (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3) {

    if (atom1 == atom3 || atom1 == atom2 || atom3 == atom2)
        return;

    atom_info a0 = __make_atom_info(atom1);
    atom_info a1 = __make_atom_info(atom2);
    atom_info a2 = __make_atom_info(atom3);

    BondAngle a_new;

    if ( a0 < a2 ) {
        a_new = make_tuple(a0, a1, a2);
    } else {
        a_new = make_tuple(a2, a1, a0);
    }

    double angle = Geom3D::angle(atom1->crd(), atom2->crd(), atom3->crd());
    ba.push_back(make_pair(a_new, angle));
}

void MoleculeBondExtractor::addDihedral (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3, const Molib::Atom* atom4) {

    if (atom4 == atom2 || atom4 == atom1 || atom3 == atom1)
        return;

    atom_info a0 = __make_atom_info(atom1);
    atom_info a1 = __make_atom_info(atom2);
    atom_info a2 = __make_atom_info(atom3);
    atom_info a3 = __make_atom_info(atom4);
    
    BondDihedral d_new;

    if ( atom1 < atom4 ) {
        d_new = make_tuple(a0, a1, a2, a3);
    } else if ( atom4 < atom1 ) {
        d_new = make_tuple(a3, a2, a1, a0);
    }

    double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                       atom3->crd(), atom4->crd());

    bd.push_back(make_pair(d_new, dihedral));
}

void MoleculeBondExtractor::addImproper (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3, const Molib::Atom* atom4) {

    if (atom4 == atom3 || atom4 == atom1 || atom3 == atom1)
        return;

    // Order should be:
    // Lowest, 2nd Lowest, Atom2, highest
    BondDihedral i_new;

    const Molib::Atom* highest = max( max( atom1, atom3), atom4);
    const Molib::Atom* lowest  = min( min( atom1, atom3), atom4);
    const Molib::Atom* middle  = max( min( atom1, atom3), min(max(atom1,atom3),atom4));

    atom_info low  = __make_atom_info(lowest);
    atom_info a2   = __make_atom_info(atom2);
    atom_info mid  = __make_atom_info(middle);
    atom_info high = __make_atom_info(highest);

    i_new = make_tuple (low, mid, a2, high);

    double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                       atom3->crd(), atom4->crd());

    bi.push_back(make_pair(i_new, dihedral));
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
        substitutions[patom1] = patom1->size() - patom1->get_num_hydrogens();
    }    

    for ( const auto& patom1 : all_my_atoms ) {
        for ( const auto& atom2 : *patom1 ) {

            const Molib::Atom& atom1 = *patom1;

            addStretch(&atom1, &atom2);

            for ( const auto& atom3 : atom2 ) {

                // Don't back track
                if ( &atom3 == &atom1 )
                    continue;

                addAngle(&atom1, &atom2, &atom3);

                for ( const auto& atom4 : atom3 ) {

                    addDihedral(&atom1, &atom2,
                                &atom3, &atom4);
                }
                for ( const auto& atom4 : atom2 ) {

                    addImproper(&atom1, &atom2,
                                &atom3, &atom4);
                }
            }
        }
    }

    return true;
}

void MoleculeBondExtractor::binStretches( const double bond_bin_size ) {
    for ( const auto& bond : bs ) {
        const size_t bin = floor(bond.second / bond_bin_size);

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
    bs.clear();
}

void MoleculeBondExtractor::binAngles( const double angle_bin_size ) {
    for ( const auto& angle : ba ) {
        const size_t bin = floor(angle.second / angle_bin_size);

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
    ba.clear();
}

void MoleculeBondExtractor::binDihedrals( const double dihedral_bin_size ) {
    for ( const auto& dihedral : bd ) {
        const int bin = floor(dihedral.second / dihedral_bin_size);

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
    bd.clear();
}

void MoleculeBondExtractor::binImpropers( const double improper_bin_size ) {
    for ( const auto& improper : bi ) {
        const int bin = floor(improper.second / improper_bin_size);

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
    bi.clear();
}

void MoleculeBondExtractor::printStretches (std::ostream& os) const {
    for ( const auto& bond : bs ) {
        os << get<0>(bond.first) << ' '
           << get<1>(bond.first) << ' ';
        os << bond.second;
        os << "\n";
    }
}

void MoleculeBondExtractor::printAngles (std::ostream& os) const {
    for ( const auto& angle : ba ) {
        os << get<0>(angle.first) << ' '
           << get<1>(angle.first) << ' '
           << get<2>(angle.first) << ' ';
        os << Geom3D::degrees(angle.second);
        os << "\n";
    }
}

void MoleculeBondExtractor::printDihedrals (std::ostream& os) const {
    for ( const auto& dihedral : bd ) {
        os << get<0>(dihedral.first) << ' '
           << get<1>(dihedral.first) << ' '
           << get<2>(dihedral.first) << ' '
           << get<3>(dihedral.first) << ' ';
        os << Geom3D::degrees(dihedral.second);
        os << "\n";
    }
}

void MoleculeBondExtractor::printImpropers (std::ostream& os) const {
    for ( const auto& improper : bi ) {
        os << get<0>(improper.first) << ' '
           << get<1>(improper.first) << ' '
           << get<2>(improper.first) << ' '
           << get<3>(improper.first) << ' ';
        os << Geom3D::degrees(improper.second);
        os << "\n";
    }
}

void MoleculeBondExtractor::printStretchBins (std::ostream& os) const {
    for ( const auto& bond : bsc ) {
        os << get<0>(bond.first.first) << " ";
        os << get<1>(bond.first.first) << " ";
        os << bond.first.second << " ";
        os << bond.second;
        os << "\n";
    }

}

void MoleculeBondExtractor::printAngleBins (std::ostream& os) const {
    for ( const auto& angle : bac ) {
        os << get<0>(angle.first.first) << " ";
        os << get<1>(angle.first.first) << " ";
        os << get<2>(angle.first.first) << " ";
        os << angle.first.second << " ";
        os << angle.second;
        os << "\n";
    }

}

void MoleculeBondExtractor::printDihedralBins (std::ostream& os) const {
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

void MoleculeBondExtractor::printImproperBins (std::ostream& os) const {
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

std::ostream& operator<< (std::ostream& os, const MoleculeBondExtractor::atom_info& ai) {
        os << help::idatm_unmask[ ai.idatm_type ] << ' ';
        os << ai.ring_size << ' ';
        os << ai.substitutions << ' ';
        return os;
}
