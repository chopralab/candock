#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/coordinate.hpp"
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "pdbreader/molecule.hpp"
#include "OpenMM.h"
using namespace std;

namespace Molib {
	class Atom;
	class Molecule;
};

namespace OMMIface {
	class ForceField;

	struct Topology {
		typedef set<pair<Molib::Atom*, Molib::Atom*>> BondedExclusions;
		typedef vector<pair<Molib::Atom*, Molib::Atom*>> Bonds;
		typedef vector<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*>> Angles;
		typedef vector<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*, Molib::Atom*>> Dihedrals;
		Molib::Atom::Vec atoms;
		Bonds bonds;
		Angles angles;
		Dihedrals dihedrals, impropers;
		BondedExclusions bonded_exclusions;
		map<const Molib::Atom*, const int> atom_to_type;

		~Topology() { dbgmsg("calling destructor of Topology"); }

		void init_topology(const Molib::Atom::Vec&, const ForceField&);
		
		
		friend ostream& operator<< (ostream& stream, const BondedExclusions &bonds);
		friend ostream& operator<< (ostream& stream, const Bonds &bonds);
		friend ostream& operator<< (ostream& stream, const Angles &angles);
		friend ostream& operator<< (ostream& stream, const Dihedrals &dihedrals);

	};

};
#endif
