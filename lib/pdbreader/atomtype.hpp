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

using namespace std;

namespace Molib {	
	class Atom;
	
	class AtomType {
	public:
		static void compute_idatm_type(Molecule &molecule);
		static void refine_idatm_type(Molecule &molecule);
		static void compute_gaff_type(Molecule &molecule);
		static void compute_ring_type(Molecule &molecule);
	};
};
#endif
