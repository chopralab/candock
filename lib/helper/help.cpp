#include "help.hpp"

namespace help {

	ostream& operator<<(ostream& os, const smiles& edges) {
		for (auto &e : edges)
			os << "vertex1 = " << e.atom_property1 
				<< " vertex2 = " << e.atom_property2
				<< " bond_property = " << e.bond_property 
				<< endl;
		return os;
	}
};
