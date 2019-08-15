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

#ifndef POSES_H
#define POSES_H

#include "molib/grid.hpp"
#include "geom3d/geom3d.hpp"
#include "segment.hpp"
#include "seed.hpp"

namespace Molib {
	class Atom;
};

namespace Linker {
	class State;
	
	class Poses {
		class AtomPoint {
			const Geom3D::Point &__crd;
			const Molib::Atom &__atom; 
			State &__state;
		public:
			AtomPoint(const Geom3D::Point &crd, const Molib::Atom &atom, State &state) : __crd(crd), __atom(atom), __state(state) {}
			const Geom3D::Point& crd() const { return __crd; }
			State& get_state() { return __state; }
			const Molib::Atom& get_atom() const { return __atom; }
			void distance(double) const {} // just dummy : needed by grid
			double radius() const { return __atom.radius(); }
			
			typedef vector<unique_ptr<AtomPoint>> UPVec;
			typedef vector<AtomPoint*> PVec;
			typedef ::Grid<AtomPoint> Grid;
		};
		
		map<Segment::Id, AtomPoint::UPVec> __atompoints;
		map<Segment::Id, AtomPoint::Grid> __grid;
		
	public:
		Poses(const Seed::Graph &seed_graph);
		State::Set get_join_states(const State &state, Segment &segment2, Molib::Atom::Pair &jatoms, 
			const double max_linker_length, const double lower_tol_seed_dist, const double upper_tol_seed_dist);
	};
}
#endif
