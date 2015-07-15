#ifndef POSES_H
#define POSES_H

#include "pdbreader/grid.hpp"

namespace Molib {
	class Atom;
};

namespace Linker {
	class State;
	
	class Poses {
		class AtomPoint {
			Molib::Atom &__atom; 
			Linker::State &__state;
		public:
			AtomPoint(Molib::Atom &atom, Linker::State &state) : __atom(atom), __state(state) {}
			Geom3D::Point& crd() { return __atom.crd(); }
			Linker::State& get_state() { return __state; }
			Molib::Atom& get_atom() { return __atom; }

			typedef vector<unique_ptr<AtomPoint>> UPVec;
			typedef vector<AtomPoint*> PVec;
		};
		
		AtomPoint::UPVec __atompoints;
		Grid<AtomPoint> __grid;
		
		Grid<AtomPoint> __initialize_grid(const State::Vec &states);
	public:
		Poses(const State::Vec &states);
		State::Set get_clashed_states(const State &state) const;
		State::Set get_join_states(const State &state) const;
	};
}
#endif
