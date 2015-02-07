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

	//~ typedef Glib::Graph<Atom> MolGraph;
	//~ typedef map<Atom*, string> AtomToRename;
	typedef vector<Atom*> AtomVec;
	typedef set<Atom*> AtomSet;
	typedef map<int, Atom*> AtomMatch;
	typedef vector<AtomMatch> AtomMatchVec;
	typedef set<Atom*> AtomSet;
	//~ typedef map<string, set<AtomSet>> Fragments;
	typedef map<string, set<AtomSet>> Fragments;
	typedef AtomSet Ring;
	typedef set<Ring> Rings;

	class Fragmenter {
	//~ public:
	//~ private:
		//~ MolGraph __mol_graph;
		AtomVec __atoms;
		const int __min_rigid_atoms, __min_gaff_group_atoms;
		void __merge_small_fragments_with_rings(AtomSet&, int);
		//~ Rings __bond_rings_to_atom_rings(const set<BondSet> &brings);
		AtomMatch __convert_to_atom_match(const map<Bond*, Bond*> &bond_match, 
			bool reverse=false);
	public:
		Fragmenter(const AtomVec &atoms);
		//~ Fragmenter(Model &model);
		//~ Fragmenter(Residue &residue);
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
		//~ void apply_rule(const Pattern &patt, const vector<string> &bond_rules);
		//~ Bonds identify_rotatable(const help::rotatable_bonds&);
		//~ Fragments identify_overlapping_rigid_segments(const Bonds&);
		//~ Fragments identify_overlapping_rigid_segments();
		//~ Fragments identify_overlapping_rigid_segments(const AtomSet &atoms);
		Fragments identify_overlapping_rigid_segments(const AtomVec &atoms);
		//~ AtomToRename rename(const help::rename_rules&);
		void flip_conjugated_gaff_types(const AtomVec &atoms);
	};
	ostream& operator<<(ostream& os, const Fragments& fragments);
	ostream& operator<<(ostream& os, const AtomMatch& atom_match);
	ostream& operator<<(ostream& os, const map<Bond*, Bond*>& bond_match);
}
#endif
