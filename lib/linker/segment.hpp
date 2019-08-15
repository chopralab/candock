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

#ifndef SEGMENT_H
#define SEGMENT_H
#include "helper/debug.hpp"
#include "molib/it.hpp"
#include "fragmenter/fragmenter.hpp"
#include "geom3d/coordinate.hpp"
#include "molib/internal.hpp"
#include "molib/molecule.hpp"
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
		typedef int Id;
		
	private:
		map<const Molib::Atom*, int> __amap;
		const Molib::Atom::Vec __atoms;
		int __seed_id; // != -1 if segment is a seed
		Id __id;
		vector<unique_ptr<State>> __state; // only seed states here!
		ConstSet __adjacent_seed_segments;
		map<const Segment*, double> __max_linker_length;
		map<const Segment*, Molib::Bond> __bond;
		map<const Segment*, Segment*> __next;
		vector<bool> __join_atom, __common_atom;

		Segment(const Molib::Atom::Vec atoms, const int &seed_id, const Segment::Id idx);
		static Paths __find_paths(const vector<unique_ptr<Segment>> &segments);
		static void __set_branching_rules(const Paths &paths);
		static bool __link_adjacent(const Graph::Path &path);
		static void __init_max_linker_length(const Paths &paths);
		static void __compute_max_linker_length(Segment::Graph::Path &path);
	public:

		int get_seed_id() const { return __seed_id; }
		bool has_next(const Segment &goal) const { return __next.count(&goal) != 0; }
		Segment &get_next(const Segment &goal) const { return *__next.at(&goal); } // get next seg in the direction of goal
		void set_next(Segment &goal, Segment &next) { __next.insert({&goal, &next}); }
		const ConstSet& get_adjacent_seed_segments() const { return __adjacent_seed_segments; }
		void set_adjacent_seed_segments(Segment &seed_seg) { __adjacent_seed_segments.insert(&seed_seg); }
		bool is_seed_adjacent(const Segment &other) const { return __adjacent_seed_segments.count(&other) != 0; }
		const Molib::Atom::Vec& get_atoms() const { return __atoms; }
		const Molib::Atom& get_atom(const int i) const { return *__atoms[i]; }
		bool has_atom(const Molib::Atom &atom) const { return __amap.count(&atom) != 0; }
                int get_idx(const Molib::Atom &atom) const { return __amap.at(&atom); }
		string get_label() const { stringstream ss; ss << *this; return ss.str(); } // graph ostream operator
		int weight() const { return 0; } // dummy for graph ostream operator
		void add_state(unique_ptr<State> s) { __state.push_back(std::move(s)); }
		bool is_seed() const { return __seed_id != -1; }
		bool is_leaf() const { return size() == 1; }
		bool is_branch() const { return size() > 2; }
		const vector<unique_ptr<State>>& get_states() const { return __state; }
		State& get_first_state() const { return *__state[0]; }
		State& get_last_state() const { return *__state.back(); }
		bool is_adjacent(const Segment &other) const { for (auto &adj : *this) if (&adj == &other) return true; return false; }
		double get_max_linker_length(const Segment &other) const { return __max_linker_length.at(&other); }
		void set_max_linker_length(const Segment &other, const double d) {	if (d > __max_linker_length[&other]) __max_linker_length[&other] = d; }
		const Molib::Bond& get_bond(const Segment &other) const { return __bond.at(&other); }
		void set_bond(const Segment &other, Molib::Atom &a1, Molib::Atom &a2);
		int adjacent_in_segment(const Molib::Atom &atom, const Molib::Atom &forbidden) const;
        Id get_id() const { return __id; }
		void set_join_atom(const Molib::Atom &atom) { __join_atom[get_idx(atom)] = true; }
		bool is_join_atom(const int i) const { return __join_atom[i]; }
		void set_common_atom(const Molib::Atom &atom) { __common_atom[get_idx(atom)] = true; }
		bool is_common_atom(const int i) const { return __common_atom[i]; }
		friend ostream& operator<< (ostream& stream, const Segment& s);
		
		static Graph create_graph(const Molib::Molecule &molecule);
	};
};
#endif
