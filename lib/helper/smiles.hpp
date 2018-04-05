#ifndef SMILES_H
#define SMILES_H

#include <vector>
#include <string>
#include <ostream>

using namespace std;

namespace help {
	struct edge { 
		string atom_property1;
		string atom_property2;
		string bond_property;
        string bond_stereo;
	};
	typedef vector<edge> smiles;
	struct rename_rule {
		smiles pattern;
		vector<string> rule;
	};
	typedef vector<rename_rule> rename_rules;

	ostream& operator<<(ostream& os, const smiles& edges);
	ostream& operator<<(ostream& os, const rename_rule& rule);
}

#endif
