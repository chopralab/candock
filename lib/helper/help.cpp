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

	ostream& operator<<(ostream& os, const rename_rule& rule) {
		os << "pattern = " << rule.pattern;
		for (auto &txt : rule.rule) 
			os << " txt = " << txt;
		os << endl;
		return os;
	}

};
