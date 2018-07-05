#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "candock/geometry/coordinate.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/help.hpp"
#include "candock/molib/molecule.hpp"
using namespace std;

namespace candock {

namespace molib {
	class Atom;
	class Molecule;
};

namespace OMMIface {
	struct ForceField;

	class Topology {
	public:
		typedef set<pair<molib::Atom*, molib::Atom*>> BondedExclusions;
		typedef vector<pair<molib::Atom*, molib::Atom*>> Bonds;
		typedef vector<tuple<molib::Atom*, molib::Atom*, molib::Atom*>> Angles;
		typedef vector<tuple<molib::Atom*, molib::Atom*, molib::Atom*, molib::Atom*>> Dihedrals;
		molib::Atom::Vec atoms;
		Bonds bonds;
		Angles angles;
		Dihedrals dihedrals, impropers;
		BondedExclusions bonded_exclusions;
	private:
		map<const molib::Atom*, const int> atom_to_type;
		map<const molib::Atom*, const int> atom_to_index;

	public:
		~Topology() { dbgmsg("calling destructor of Topology"); }

		Topology& add_topology(const molib::Atom::Vec &atoms, const ForceField &ffield);
		int get_index(const molib::Atom &atom) const;
		int get_type(const molib::Atom &atom) const;
		
	};

	ostream& operator<< (ostream& stream, const Topology::BondedExclusions &bonds);
	ostream& operator<< (ostream& stream, const Topology::Bonds &bonds);
	ostream& operator<< (ostream& stream, const Topology::Angles &angles);
	ostream& operator<< (ostream& stream, const Topology::Dihedrals &dihedrals);
};

}

#endif
