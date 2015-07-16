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
			Molib::Atom &__atom; 
			State &__state;
		public:
			AtomPoint(Molib::Atom &atom, State &state) : __atom(atom), __state(state) {}
			Geom3D::Point& crd() { return __atom.crd(); }
			State& get_state() { return __state; }
			Molib::Atom& get_atom() { return __atom; }
			void distance(double d) const {} // just dummy : needed by grid

			typedef vector<unique_ptr<AtomPoint>> UPVec;
			typedef vector<AtomPoint*> PVec;
		};
		
		//~ AtomPoint::UPVec __atompoints, __joinatompoints;
		//~ map<Segment::Id, Grid<AtomPoint>> __grid, __grid_join_atoms;
		//~ 
		AtomPoint::UPVec __atompoints;
		map<Segment::Id, Grid<AtomPoint>> __grid;
		const double __tol_seed_dist;
		
		void __initialize_grid(const Seed::Graph &seed_graph);
	public:
		Poses(const Seed::Graph &seed_graph, const double tol_seed_dist);
		State::Set get_clashed_states(const State &state, Segment &segment2);
		State::Set get_join_states(const State &state, Segment &segment2, Molib::Atom::Pair &jatoms, const double max_linker_length);
	};
}
#endif
