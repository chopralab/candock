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

#include "smiles.hpp"

namespace help {
	ostream& operator<<(ostream& os, const smiles& edges) {
		for (auto &e : edges)
			os << '{' << e.atom_property1 
				<< ',' << e.atom_property2
				<< ',' << e.bond_property << '}';
		return os;
	}

	ostream& operator<<(ostream& os, const rename_rule& rule) {
		os << "PATTERN = " << rule.pattern << " RULE = {";
		for (auto &txt : rule.rule) os << txt << ',';
		os << '}';
		//~ os << endl;
		return os;
	}
}
