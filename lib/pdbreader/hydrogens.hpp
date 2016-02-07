#ifndef HYDROGENS_H
#define HYDROGENS_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <cstdlib>
#include "it.hpp"
#include "atom.hpp"

using namespace std;

namespace Molib {	
	class Atom;
	class Hydrogens {
	public:
		static void compute_hydrogen(const Atom::Vec &atoms);
		static void erase_hydrogen(const Atom::Vec &atoms);
	};
};
#endif
