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

namespace Molib {
	class State;
	class Segment;
	typedef set<Segment*> SegmentSet;
	typedef set<const Segment*> ConstSegmentSet;
	typedef vector<Segment*> SegmentVec;
	typedef pair<const Segment*, const Segment*> ConstSegPair;
	class Segment : public template_vector_container<Segment*, Segment> {
		const AtomSet &__atoms;
		const AtomSet *__seed_atoms;
		string __name; // name of segment is seed name or empty string
		vector<unique_ptr<State>> __state; // only seed states here!
		ConstSegmentSet __adjacent_seed_segments;
		map<const Segment*, const double> __max_linker_length;
		map<const Segment*, Bond> __bond;
		map<const Segment*, Segment*> __next;
	public:
		Segment(const AtomSet &atoms) : __atoms(atoms), __name("") {}
		string get_name() const { return __name; }
		bool has_next(const Segment &goal) const { return __next.count(&goal); }
		Segment &get_next(const Segment &goal) const { return *__next.at(&goal); } // get next seg in the direction of goal
		void set_next(Segment &goal, Segment &next) { __next.insert({&goal, &next}); }
		const AtomSet& get_seed_atoms() const { return *__seed_atoms; }
		const ConstSegmentSet& get_adjacent_seed_segments() const { return __adjacent_seed_segments; };
		void set_adjacent_seed_segments(Segment &seed_seg) { __adjacent_seed_segments.insert(&seed_seg); };
		bool is_seed_adjacent(const Segment &other) const { return __adjacent_seed_segments.count(&other); }
		const AtomSet& get_atoms() const { return __atoms; }
		bool has_atom(Atom &atom) const { return __atoms.count(&atom); }
		const string get_label() const { stringstream ss; ss << *this; return ss.str(); } // graph ostream operator
		const int weight() const { return 0; } // dummy for graph ostream operator
		void set_name(string nm) { __name = nm; }
		void set_seed_atoms(const AtomSet &sa) { __seed_atoms = &sa; }
		void add_state(unique_ptr<State> s) { __state.push_back(std::move(s)); }
		bool is_seed() const { return !__name.empty(); }
		bool is_leaf() const { return size() == 1; }
		bool is_branch() const { return size() > 2; }
		const vector<unique_ptr<State>>& get_states() const { return __state; }
		State& get_first_state() const { return *__state[0]; }
		bool is_adjacent(const Segment &other) const { for (auto &adj : *this) if (&adj == &other) return true; return false; }
		double get_max_linker_length(const Segment &other) const { return __max_linker_length.at(&other); }
		void set_max_linker_length(const Segment &other, const double d) {	__max_linker_length.insert({&other, d}); }
		const Bond& get_bond(const Segment &other) const { return __bond.at(&other); }
		void set_bond(const Segment &other, const Bond &b) { __bond.insert({&other, b}); }
		//~ const Atom& adjacent_in_segment(const Atom &atom, const Atom &forbidden) const { for (auto &adj : atom) if (&adj != &forbidden && has_atom(adj)) return adj; }
		const Atom& adjacent_in_segment(const Atom &atom, const Atom &forbidden) const;
		friend ostream& operator<< (ostream& stream, const Segment& s) {
			stream << "Segment(" << s.__name << ") : atom numbers = ";
			for (auto &pa : s.__atoms) stream << pa->atom_number() << " ";
			return stream;
		}
	};
	class Seed : public template_vector_container<Seed*, Seed> {
		Segment &__seg;
	public:
		Seed(Segment &seg) : __seg(seg) {}
		Segment& get_segment() const { return __seg; } 
		const string get_label() const { return __seg.get_label(); } // graph ostream operator
		const int weight() const { return 0; } // dummy for graph ostream operator
		friend ostream& operator<< (ostream& stream, const Seed& s) {
			return stream << "Seed = " << s.get_segment().get_name() << endl;
		}
	};
};
#endif
