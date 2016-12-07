#ifndef HYDROGENS_H
#define HYDROGENS_H

// TODO: Possible place these functions in residue.cpp
namespace Molib {
	class Atom;
	class Residue;
	class Hydrogens {
	public:
		// These functions modify the reside structure of the atoms given to them.
		// thus I feel it is correct to pass this into the functions, not the atom
		// vectors.
		static void compute_hydrogen(Residue &res);
		static void erase_hydrogen(Residue &res);
	};
};
#endif
