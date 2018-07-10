#ifndef FRAGMENTER_H
#define FRAGMENTER_H
#include "candock/helper/smiles.hpp"
#include "candock/molib/bond.hpp"
#include "candock/molib/atom.hpp"

#include <map>
#include <vector>
#include <set>

namespace candock {

namespace molib {
	class Unique;

	typedef std::map<int, Atom*> AtomMatch;
	typedef std::vector<AtomMatch> AtomMatchVec;
	typedef Atom::Set Ring;
	typedef std::set<Ring> Rings;

	class Fragmenter {
	public:
	
		class Fragment {
			Atom::Set __core, __join;
			int __seed_id;
		public:
			Fragment(Atom::Set core, Atom::Set join, Unique &u);
			Fragment(Atom::Set core, Atom::Set join, int seed_id) : __core(core), __join(join), __seed_id(seed_id) {}
			Atom::Set& get_core() { return __core; }
			Atom::Set& get_join() { return __join; }
			const Atom::Set& get_core() const { return __core; }
			const Atom::Set& get_join() const { return __join; }
			Atom::Set get_all() const { Atom::Set all(__core); all.insert(__join.begin(), __join.end()); return all; }
			bool is_seed() const { return __seed_id != -1; }
			int get_seed_id() const { return __seed_id; }
			int size() const { return __core.size() + __join.size(); }

			typedef std::vector<Fragment> Vec;
			friend std::ostream& operator<<(std::ostream& os, const Vec& fragments);
		};
	
	private:
		Atom::Vec __atoms;
		AtomMatch __convert_to_atom_match(const std::map<Bond*, Bond*> &bond_match, 
			bool reverse=false);
	public:
		Fragmenter(const Atom::Vec &atoms);
		Rings identify_rings();
		Rings identify_fused_rings();
		AtomMatchVec grep(const help::smiles &smi);
		void apply_rule(const AtomMatch &m, 
			const std::vector<std::string> &rules, Atom::Set &visited); // atom rule
		void apply_rule(const AtomMatch &m, 
			const std::vector<std::string> &rules, BondSet &visited); // bond rule
		void substitute_bonds(const help::rename_rules &rrules);
		void substitute_atoms(const help::rename_rules &rrules);
		Fragment::Vec identify_overlapping_rigid_segments(const Atom::Vec &atoms, Unique &u);
		void flip_conjugated_gaff_types(const Atom::Vec &atoms);
	};
	std::ostream& operator<<(std::ostream& os, const AtomMatchVec& atom_match_vec);
	std::ostream& operator<<(std::ostream& os, const AtomMatch& atom_match);
	std::ostream& operator<<(std::ostream& os, const std::map<Bond*, Bond*>& bond_match);
}

}

#endif
