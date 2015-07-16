#include "poses.hpp"
#include "state.hpp"
#include "segment.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "graph/mcqd.hpp"
#include <queue>

using namespace std;

namespace Linker {
	Poses::Poses(const Segment::Vec &segments) {
		try {
			__initialize_grid(segments);
		} catch (...) {
			dbgmsg("FAILURE: constructor of Poses failed...");
			throw;
		}
	}
	
	void Poses::__initialize_grid(const Segment::Vec &segments) {
		for (auto &segment : segments) {
			for (auto &state : segment.get_states()) {
				for (auto &atom : state.get_atoms()) {
					__atompoints[segment.get_id].push_back(unique_ptr<AtomPoint>(new AtomPoint(atom, state)));
				}
				for (auto &atom : state.get_join_atoms()) {
					__joinatompoints[segment.get_id()].push_back(unique_ptr<AtomPoint>(new AtomPoint(atom, state)));
				}
			}
			__grid[segment.get_id()] = Grid<AtomPoint>(__atompoints[segment.get_id()]);
			__grid_join_atoms[segment.get_id()] = Grid<AtomPoint>(__joinatompoints[segment.get_id()]);
		}
	}
	
	State::Set Poses::get_clashed_states(const State &state, const Segment &segment2) const {
		State::Set clashed;
		Grid<AtomPoint> &g = __grid[segment2.get_id()];
		for (auto &atom : state.get_atoms()) {
			AtomPoint::PVec neighbors = g.get_neighbors(atom.crd(), atom.radius());
			for (auto &patompoint : neighbors) {
				clashed.insert(&patompoint->get_state());
		}
		return clashed;
	}

	State::Set get_join_states(const State &state, const Segment &segment2, const Segment::JoinAtoms &jatoms, const double max_linker_length) const {
		State::Set join;
		Grid<AtomPoint> &g = __grid_join_atoms[segment2];
		for (auto &atom : state.get_join_atoms()) {
			AtomPoint::PVec neighbors = g.get_neighbors_within_tolerance(atom.crd(), 
				max_linker_length, tolerance);
					
			for (auto &patompoint : neighbors) {
				if (jatoms.first == atom && jatoms.second == patompoint->get_atom()) {
					join.insert(&patompoint->get_state());
				}
			}
		}
		return join;
	}
};
