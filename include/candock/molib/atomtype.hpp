#ifndef ATOMTYPE_H
#define ATOMTYPE_H
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "candock/molib/atom.hpp"
#include "candock/molib/it.hpp"

namespace candock {

namespace molib {
class Atom;
class Molecule;

namespace AtomType {
void compute_idatm_type(const Atom::Vec& atoms);
void refine_idatm_type(const Atom::Vec& atoms);
void compute_gaff_type(const Atom::Vec& atoms);
void compute_ring_type(const Atom::Vec& atoms);

std::tuple<double, size_t, size_t, size_t> determine_lipinski(
    const Atom::Vec& atoms);
};
};
}

#endif
