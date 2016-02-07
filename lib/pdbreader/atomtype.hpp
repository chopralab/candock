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
	
	class AtomType {
	public:
		static void compute_idatm_type(const Atom::Vec &atoms);
		static void refine_idatm_type(const Atom::Vec &atoms);
		static void compute_gaff_type(const Atom::Vec &atoms);
		static void compute_ring_type(const Atom::Vec &atoms);
	};
};
#endif
