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
		int get_index(const Molib::Atom &atom) const { return atom_to_index.at(&atom); }
		int get_type(const Molib::Atom &atom) const { return atom_to_type.at(&atom); }
		
	};

	ostream& operator<< (ostream& stream, const Topology::BondedExclusions &bonds);
	ostream& operator<< (ostream& stream, const Topology::Bonds &bonds);
	ostream& operator<< (ostream& stream, const Topology::Angles &angles);
	ostream& operator<< (ostream& stream, const Topology::Dihedrals &dihedrals);


};
#endif
