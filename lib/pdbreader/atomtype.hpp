#ifndef ATOMTYPE_H
#define ATOMTYPE_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "it.hpp"
#include "atom.hpp"

using namespace std;

namespace Molib {
	class Atom;
	class Molecule;
	
	namespace AtomType {
		void compute_idatm_type(const Atom::Vec &atoms);
		void refine_idatm_type(const Atom::Vec &atoms);
		void compute_gaff_type(const Atom::Vec &atoms);
		void compute_ring_type(const Atom::Vec &atoms);
                
                std::tuple<double, int, int> determine_lipinski(const Atom::Vec &atoms);
	};
};
#endif
