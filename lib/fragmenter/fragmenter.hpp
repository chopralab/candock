#ifndef FRAGMENTER_H
#define FRAGMENTER_H
#include "helper/debug.hpp"
#include "graph/graph.hpp"
#include "helper/help.hpp"
#include "pdbreader/bond.hpp"
#include <tuple>
#include <vector>
#include <map>
#include <set>
#include <functional>
using namespace std;
namespace Molib {
	class Unique;
	class Atom;
	class Residue;
	class Model;
	class Molecule;

	typedef vector<Atom*> AtomVec;
	typedef set<Atom*> AtomSet;
	typedef map<int, Atom*> AtomMatch;
	typedef vector<AtomMatch> AtomMatchVec;
	typedef set<Atom*> AtomSet;
	typedef map<string, set<AtomSet>> Fragments;
	typedef AtomSet Ring;
	typedef set<Ring> Rings;

	class Fragmenter {
		AtomVec __atoms;
		const int __min_rigid_atoms, __min_gaff_group_atoms;
		void __merge_small_fragments_with_rings(AtomSet&, int);
		AtomMatch __convert_to_atom_match(const map<Bond*, Bond*> &bond_match, 
			bool reverse=false);
	public:
		Fragmenter(const AtomVec &atoms);
		Rings identify_rings();
		Rings identify_fused_rings();
		Fragments identify_seeds(const Fragments &rigid, Unique &u);
		AtomMatchVec grep(const help::smiles &smi);
		void apply_rule(const AtomMatch &m, 
			const vector<string> &rules, AtomSet &visited); // atom rule
		void apply_rule(const AtomMatch &m, 
			const vector<string> &rules, BondSet &visited); // bond rule
		void substitute_bonds(const help::rename_rules &rrules);
		void substitute_atoms(const help::rename_rules &rrules);
		Fragments identify_overlapping_rigid_segments(const AtomVec &atoms);
		void flip_conjugated_gaff_types(const AtomVec &atoms);
	};
	ostream& operator<<(ostream& os, const Fragments& fragments);
	ostream& operator<<(ostream& os, const AtomMatchVec& atom_match_vec);
	ostream& operator<<(ostream& os, const AtomMatch& atom_match);
	ostream& operator<<(ostream& os, const map<Bond*, Bond*>& bond_match);
}
#endif
