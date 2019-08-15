/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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
