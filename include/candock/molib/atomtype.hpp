#ifndef ATOMTYPE_H
#define ATOMTYPE_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "candock/molib/it.hpp"
#include "candock/molib/atom.hpp"

namespace Molib {
	class Atom;
	class Molecule;
	
	namespace AtomType {
		void compute_idatm_type(const Atom::Vec &atoms);
		void refine_idatm_type(const Atom::Vec &atoms);
		void compute_gaff_type(const Atom::Vec &atoms);
		void compute_ring_type(const Atom::Vec &atoms);
                
                std::tuple<double, size_t, size_t, size_t> determine_lipinski(const Atom::Vec &atoms);
	};
};
#endif
