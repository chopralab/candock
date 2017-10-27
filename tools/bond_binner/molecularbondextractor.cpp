#include "molecularbondextractor.hpp"

namespace AtomInfo {

int MolecularBondExtractor::__ring_size( const Molib::Atom* atom ) const {
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

void MolecularBondExtractor::__print_atom( const Molib::Atom* atom, std::ostream& os, char delim ) const {
    os << help::idatm_unmask[ atom->idatm_type() ] << delim;
    os << __ring_size(atom) << delim;
    os << substitutions.at(atom) << delim;
}

atom_info MolecularBondExtractor::__make_atom_info(const Molib::Atom* atom) {
    atom_info a = { atom->idatm_type(),
                    __ring_size(atom),
                    substitutions.at(atom)
    };

    return a;
}

void MolecularBondExtractor::reset() {
    substitutions.clear();
    all_rings.clear();
    all_rings_sizes.clear();
}

void MolecularBondExtractor::addStretch (const Molib::Atom* atom1, const Molib::Atom* atom2) {

    // Avoid a double count if starting with a later atom
    if (atom1 >= atom2)
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

void MolecularBondExtractor::addAngle (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3) {

    // Middle and end cases not really possible, but the check is cheap
    if (atom1 >= atom3 || atom1 == atom2 || atom3 == atom2)
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

void MolecularBondExtractor::addDihedral (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3, const Molib::Atom* atom4) {

    // Some later checks are probably not needed
    if ( atom1 >= atom4 || atom3 == atom1 || atom4 == atom2 || atom2 == atom3)
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

void MolecularBondExtractor::addImproper (const Molib::Atom* atom1, const Molib::Atom* atom2, const Molib::Atom* atom3, const Molib::Atom* atom4) {

    // Leave atom2 alone, it is the pivot atom
    if (atom1 >= atom3 || atom1 >= atom4 || atom3 >= atom4)
        return;

    // Order should be:
    // Lowest, 2nd Lowest, Atom2, highest
    BondDihedral i_new;

    atom_info a0 = __make_atom_info(atom1);
    atom_info a1 = __make_atom_info(atom2);
    atom_info a2 = __make_atom_info(atom3);
    atom_info a3 = __make_atom_info(atom4);

    atom_info lowest  = min( min( a0, a1), a3);
    atom_info middle  = max( min( a0, a2), min(max(a0,a2),a3));
    atom_info highest = max( max( a0, a2), a3);

    i_new = make_tuple (lowest, middle, a1, highest);

    double dihedral = Geom3D::dihedral(atom1->crd(), atom2->crd(),
                                       atom3->crd(), atom4->crd());

    bi.push_back(make_pair(i_new, dihedral));
}


bool MolecularBondExtractor::addMolecule( const Molib::Molecule& mol ) {

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

const VBondStretches& MolecularBondExtractor::stretches() const {
    return bs;
}

const VBondAngles&    MolecularBondExtractor::angles() const {
    return ba;
}
const VBondDihedrals& MolecularBondExtractor::dihedrals() const {
    return bd;
}

const VBondDihedrals& MolecularBondExtractor::impropers() const {
    return bi;
}

void MolecularBondExtractor::printStretches (std::ostream& os) const {
    for ( const auto& bond : bs ) {
        os << get<0>(bond.first) << ' '
           << get<1>(bond.first) << ' ';
        os << bond.second;
        os << "\n";
    }
}

void MolecularBondExtractor::printAngles (std::ostream& os) const {
    for ( const auto& angle : ba ) {
        os << get<0>(angle.first) << ' '
           << get<1>(angle.first) << ' '
           << get<2>(angle.first) << ' ';
        os << Geom3D::degrees(angle.second);
        os << "\n";
    }
}

void MolecularBondExtractor::printDihedrals (std::ostream& os) const {
    for ( const auto& dihedral : bd ) {
        os << get<0>(dihedral.first) << ' '
           << get<1>(dihedral.first) << ' '
           << get<2>(dihedral.first) << ' '
           << get<3>(dihedral.first) << ' ';
        os << Geom3D::degrees(dihedral.second);
        os << "\n";
    }
}

void MolecularBondExtractor::printImpropers (std::ostream& os) const {
    for ( const auto& improper : bi ) {
        os << get<0>(improper.first) << ' '
           << get<1>(improper.first) << ' '
           << get<2>(improper.first) << ' '
           << get<3>(improper.first) << ' ';
        os << Geom3D::degrees(improper.second);
        os << "\n";
    }
}

std::ostream& operator<< (std::ostream& os, const atom_info& ai) {
        os << help::idatm_unmask[ ai.idatm_type ] << ' ';
        os << ai.ring_size << ' ';
        os << ai.substitutions << ' ';
        return os;
}

}
