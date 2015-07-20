#ifndef SEGMENT_H
#define SEGMENT_H
#include "helper/debug.hpp"
#include "pdbreader/it.hpp"
#include "fragmenter/fragmenter.hpp"
#include "geom3d/coordinate.hpp"
#include "pdbreader/internal.hpp"
#include "pdbreader/molecule.hpp"
#include "graph/graph.hpp"
#include <tuple>
#include <functional>
#include "state.hpp"

namespace Linker {

	class Segment : public template_vector_container<Segment*, Segment> {
	public:

		typedef set<Segment*> Set;
		typedef set<const Segment*> ConstSet;
		typedef vector<Segment*> Vec;
		typedef pair<const Segment*, const Segment*> ConstPair;
		typedef Glib::Graph<Segment> Graph;
		typedef map<ConstPair, Graph::Path> Paths;
		typedef Segment* Id;
		
	private:
		map<const Molib::Atom*, int> __amap;
		const Molib::Atom::Vec __atoms;
		int __seed_id; // != -1 if segment is a seed
		vector<unique_ptr<State>> __state; // only seed states here!
		ConstSet __adjacent_seed_segments;
		map<const Segment*, const double> __max_linker_length;
		map<const Segment*, Molib::Bond> __bond;
		map<const Segment*, Segment*> __next;
		
	public:

		//~ Segment(const Molib::Atom::Set atoms, const int &seed_id) : __atoms(atoms), __seed_id(seed_id) {}
		//~ Segment(const Molib::Atom::Set xxxatoms, const int &seed_id) : __xxxatoms(xxxatoms), __seed_id(seed_id) { __atoms.assign(xxxatoms.begin(), xxxatoms.end()); }
		Segment(const Molib::Atom::Vec atoms, const int &seed_id) : __atoms(atoms), __seed_id(seed_id) { for (int i = 0; i < atoms.size(); ++i) __amap[atoms[i]] = i; }
		int get_seed_id() const { return __seed_id; }
		bool has_next(const Segment &goal) const { return __next.count(&goal); }
		Segment &get_next(const Segment &goal) const { return *__next.at(&goal); } // get next seg in the direction of goal
		void set_next(Segment &goal, Segment &next) { __next.insert({&goal, &next}); }
		const ConstSet& get_adjacent_seed_segments() const { return __adjacent_seed_segments; };
		void set_adjacent_seed_segments(Segment &seed_seg) { __adjacent_seed_segments.insert(&seed_seg); };
		bool is_seed_adjacent(const Segment &other) const { return __adjacent_seed_segments.count(&other); }
		//~ const Molib::Atom::Set& get_atoms() const { return __atoms; }
		const Molib::Atom::Vec& get_atoms() const { return __atoms; }
		const Molib::Atom& get_atom(const int i) const { return *__atoms.at(i); }
		bool has_atom(const Molib::Atom &atom) const { return __amap.count(&atom); }
		//~ bool has_atom(Molib::Atom &atom) const { return get_idx(atom) <__atoms.size(); }
		//~ bool has_atom(const int i) const { return i < __atoms.size(); }
		//~ int get_idx(const Molib::Atom &atom) const { for (int i = 0; i < __atoms.size(); ++i) if (&atom == __atoms[i]) return i; throw Error("die : atom not found in segment"); }
		//~ const int get_idx(const Molib::Atom &atom) const { try { return __amap.at(&atom); } catch(const out_of_range&) { throw Error("die : atom not found in segment"); } }
		const int get_idx(const Molib::Atom &atom) const { return __amap.at(&atom); }
		const string get_label() const { stringstream ss; ss << *this; return ss.str(); } // graph ostream operator
		const int weight() const { return 0; } // dummy for graph ostream operator
		void add_state(unique_ptr<State> s) { __state.push_back(std::move(s)); }
		bool is_seed() const { return __seed_id != -1; }
		bool is_leaf() const { return size() == 1; }
		bool is_branch() const { return size() > 2; }
		const vector<unique_ptr<State>>& get_states() const { return __state; }
		State& get_first_state() const { return *__state[0]; }
		State& get_last_state() const { return *__state.back(); }
		bool is_adjacent(const Segment &other) const { for (auto &adj : *this) if (&adj == &other) return true; return false; }
		double get_max_linker_length(const Segment &other) const { return __max_linker_length.at(&other); }
		void set_max_linker_length(const Segment &other, const double d) {	__max_linker_length.insert({&other, d}); }
		const Molib::Bond& get_bond(const Segment &other) const { return __bond.at(&other); }
		//~ void set_bond(const Segment &other, const Molib::Bond &b) { __bond.insert({&other, b}); }
		void set_bond(const Segment &other, Molib::Atom &a1, Molib::Atom &a2);
		//~ const Molib::Atom& adjacent_in_segment(const Molib::Atom &atom, const Molib::Atom &forbidden) const;
		const int adjacent_in_segment(const Molib::Atom &atom, const Molib::Atom &forbidden) const;
		Id get_id() { return this; }
		friend ostream& operator<< (ostream& stream, const Segment& s);
		
		static Graph create_graph(const Molib::Molecule &molecule);
	};
};
#endif
