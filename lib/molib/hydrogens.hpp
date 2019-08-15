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
