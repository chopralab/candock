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

#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/coordinate.hpp"
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "molib/molecule.hpp"
using namespace std;

namespace Molib {
	class Atom;
	class Molecule;
};

namespace OMMIface {
	struct ForceField;

	class Topology {
	public:
		typedef set<pair<Molib::Atom*, Molib::Atom*>> BondedExclusions;
		typedef vector<pair<Molib::Atom*, Molib::Atom*>> Bonds;
		typedef vector<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*>> Angles;
		typedef vector<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*, Molib::Atom*>> Dihedrals;
		Molib::Atom::Vec atoms;
		Bonds bonds;
		Angles angles;
		Dihedrals dihedrals, impropers;
		BondedExclusions bonded_exclusions;
	private:
		map<const Molib::Atom*, const int> atom_to_type;
		map<const Molib::Atom*, const int> atom_to_index;

	public:
		~Topology() { dbgmsg("calling destructor of Topology"); }

		Topology& add_topology(const Molib::Atom::Vec &atoms, const ForceField &ffield);
		int get_index(const Molib::Atom &atom) const;
		int get_type(const Molib::Atom &atom) const;
		
	};

	ostream& operator<< (ostream& stream, const Topology::BondedExclusions &bonds);
	ostream& operator<< (ostream& stream, const Topology::Bonds &bonds);
	ostream& operator<< (ostream& stream, const Topology::Angles &angles);
	ostream& operator<< (ostream& stream, const Topology::Dihedrals &dihedrals);


};
#endif
