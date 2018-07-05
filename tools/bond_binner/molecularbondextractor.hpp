#ifndef MOLECULARBONDEXTRACTOR_HPP
#define MOLECULARBONDEXTRACTOR_HPP

#include <tuple>
#include <map>
#include "candock/molib/atom.hpp"
#include "candock/molib/molecule.hpp"
#include "atominfotypes.hpp"

namespace AtomInfo {

class MolecularBondExtractor {

private:
    VBondStretchValues bs;
    VBondAngleValues    ba;
    VBondDihedralValues bd;
    VBondDihedralValues bi; // Impropers

    candock::Glib::Graph<candock::Molib::Atom>::Cycles all_rings;
    candock::Glib::Graph<candock::Molib::Atom>::VertexRingMap all_rings_sizes;
    std::map<const candock::Molib::Atom*, size_t> substitutions;

    size_t __ring_size( const candock::Molib::Atom* atom ) const;

    void __print_atom( const candock::Molib::Atom* atom, std::ostream& os, char delim ) const;

    atom_info __make_atom_info( const candock::Molib::Atom* atom);
    
    void addStretch (const candock::Molib::Atom* atom1, const candock::Molib::Atom* atom2 );
    void addAngle   (const candock::Molib::Atom* atom1, const candock::Molib::Atom* atom2, const candock::Molib::Atom* atom3);
    void addDihedral(const candock::Molib::Atom* atom1, const candock::Molib::Atom* atom2, const candock::Molib::Atom* atom3, const candock::Molib::Atom* atom4);
    void addImproper(const candock::Molib::Atom* atom1, const candock::Molib::Atom* atom2, const candock::Molib::Atom* atom3, const candock::Molib::Atom* atom4);

public:

    void reset();

    bool addMolecule( const candock::Molib::Molecule& mol );

    const VBondStretchValues& stretches() const;
    const VBondAngleValues&    angles() const;
    const VBondDihedralValues& dihedrals() const;
    const VBondDihedralValues& impropers() const;
    
    void printStretches (std::ostream& os) const;
    void printAngles (std::ostream& os) const;
    void printDihedrals (std::ostream& os) const;
    void printImpropers (std::ostream& os) const;

};

}

#endif
