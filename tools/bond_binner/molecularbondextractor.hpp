#ifndef MOLECULARBONDEXTRACTOR_HPP
#define MOLECULARBONDEXTRACTOR_HPP

#include <map>
#include <tuple>
#include "atominfotypes.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/molecule.hpp"

namespace AtomInfo {

class MolecularBondExtractor {
   private:
    VBondStretchValues bs;
    VBondAngleValues ba;
    VBondDihedralValues bd;
    VBondDihedralValues bi;  // Impropers

    candock::graph::Graph<candock::molib::Atom>::Cycles all_rings;
    candock::graph::Graph<candock::molib::Atom>::VertexRingMap all_rings_sizes;
    std::map<const candock::molib::Atom*, size_t> substitutions;

    size_t __ring_size(const candock::molib::Atom* atom) const;

    void __print_atom(const candock::molib::Atom* atom, std::ostream& os,
                      char delim) const;

    atom_info __make_atom_info(const candock::molib::Atom* atom);

    void addStretch(const candock::molib::Atom* atom1,
                    const candock::molib::Atom* atom2);
    void addAngle(const candock::molib::Atom* atom1,
                  const candock::molib::Atom* atom2,
                  const candock::molib::Atom* atom3);
    void addDihedral(const candock::molib::Atom* atom1,
                     const candock::molib::Atom* atom2,
                     const candock::molib::Atom* atom3,
                     const candock::molib::Atom* atom4);
    void addImproper(const candock::molib::Atom* atom1,
                     const candock::molib::Atom* atom2,
                     const candock::molib::Atom* atom3,
                     const candock::molib::Atom* atom4);

   public:
    void reset();

    bool addMolecule(const candock::molib::Molecule& mol);

    const VBondStretchValues& stretches() const;
    const VBondAngleValues& angles() const;
    const VBondDihedralValues& dihedrals() const;
    const VBondDihedralValues& impropers() const;

    void printStretches(std::ostream& os) const;
    void printAngles(std::ostream& os) const;
    void printDihedrals(std::ostream& os) const;
    void printImpropers(std::ostream& os) const;
};
}

#endif
