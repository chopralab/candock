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
	typedef AtomSet Ring;
	typedef set<Ring> Rings;

	class Fragmenter {
	public:
	
		class Fragment {
			AtomSet __core, __join;
			int __seed_id;
		public:
			Fragment(AtomSet core, AtomSet join, Unique &u);
			Fragment(AtomSet core, AtomSet join, int seed_id) : __core(core), __join(join), __seed_id(seed_id) {}
			AtomSet& get_core() { return __core; }
			AtomSet& get_join() { return __join; }
			const AtomSet& get_core() const { return __core; }
			const AtomSet& get_join() const { return __join; }
			AtomSet get_all() const { AtomSet all(__core); all.insert(__join.begin(), __join.end()); return all; }
			bool is_seed() const { return __seed_id != -1; }
			int get_seed_id() const { return __seed_id; }
			int size() const { return __core.size() + __join.size(); }

			typedef vector<Fragment> Vec;
			friend ostream& operator<<(ostream& os, const Vec& fragments);
		};
	
	private:
		AtomVec __atoms;
		AtomMatch __convert_to_atom_match(const map<Bond*, Bond*> &bond_match, 
			bool reverse=false);
	public:
		Fragmenter(const AtomVec &atoms);
		Rings identify_rings();
		Rings identify_fused_rings();
		AtomMatchVec grep(const help::smiles &smi);
		void apply_rule(const AtomMatch &m, 
			const vector<string> &rules, AtomSet &visited); // atom rule
		void apply_rule(const AtomMatch &m, 
			const vector<string> &rules, BondSet &visited); // bond rule
		void substitute_bonds(const help::rename_rules &rrules);
		void substitute_atoms(const help::rename_rules &rrules);
		Fragment::Vec identify_overlapping_rigid_segments(const AtomVec &atoms, Unique &u);
		void flip_conjugated_gaff_types(const AtomVec &atoms);
	};
	ostream& operator<<(ostream& os, const AtomMatchVec& atom_match_vec);
	ostream& operator<<(ostream& os, const AtomMatch& atom_match);
	ostream& operator<<(ostream& os, const map<Bond*, Bond*>& bond_match);
}
#endif
