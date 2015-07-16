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
		
		map<Segment::Id, AtomPoint::UPVec> __atompoints, __joinatompoints;
		map<Segment::Id, Grid<AtomPoint>> __grid, __grid_join_atoms;
		
		void __initialize_grid(const Segment::Vec &segments);
	public:
		Poses(const Segment::Vec &segments);
		State::Set get_clashed_states(const State &state, const Segment &segment2) const;
		State::Set get_join_states(const State &state, const Segment &segment2, const Segment::JoinAtoms &jatoms, const double max_linker_length) const;
	};
}
#endif
