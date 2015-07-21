#ifndef POSES_H
#define POSES_H

#include "pdbreader/grid.hpp"
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
			void distance(double d) const {} // just dummy : needed by grid
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
			const double max_linker_length, const double tol_seed_dist);
	};
}
#endif
