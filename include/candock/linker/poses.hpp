#ifndef POSES_H
#define POSES_H

#include "candock/molib/grid.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/linker/segment.hpp"
#include "candock/linker/seed.hpp"

namespace candock {

namespace molib {
	class Atom;
};

namespace linker {
	class State;
	
	class Poses {
		class AtomPoint {
			const geometry::Point &__crd;
			const molib::Atom &__atom; 
			State &__state;
		public:
			AtomPoint(const geometry::Point &crd, const molib::Atom &atom, State &state) : __crd(crd), __atom(atom), __state(state) {}
			const geometry::Point& crd() const { return __crd; }
			State& get_state() { return __state; }
			const molib::Atom& get_atom() const { return __atom; }
			void distance(double) const {} // just dummy : needed by grid
			double radius() const { return __atom.radius(); }
			
			typedef vector<unique_ptr<AtomPoint>> UPVec;
			typedef vector<AtomPoint*> PVec;
			typedef candock::molib::Grid<AtomPoint> Grid;
		};
		
		map<Segment::Id, AtomPoint::UPVec> __atompoints;
		map<Segment::Id, AtomPoint::Grid> __grid;
		
	public:
		Poses(const Seed::Graph &seed_graph);
		State::Set get_join_states(const State &state, Segment &segment2, molib::Atom::Pair &jatoms, 
			const double max_linker_length, const double lower_tol_seed_dist, const double upper_tol_seed_dist);
	};
}

}

#endif
