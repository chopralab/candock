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
	Poses::Poses(const State::Vec &states) {
		try {
			__grid = __initialize_grid(states);
		} catch (...) {
			dbgmsg("FAILURE: constructor of Poses failed...");
			throw;
		}
	}
	
	Grid<Poses::AtomPoint> Poses::__initialize_grid(const State::Vec &states) {
		for (auto &state : states) {
			for (state.get_atoms()) {
				__atompoints.push_back(unique_ptr<AtomPoint>(new AtomPoint(atom, state)));
			}
		}
		return Grid<AtomPoint> grid(atompoints);
	}
	
	State::Set get_join_states(const State &state, const double max_linker_length) const {
		State::Set join;
		for (auto &atom : state.get_join_atoms()) {
			for (auto &segment2 : other_segments) {
				AtomPoint::PVec neighbors = 
					__grid[segment2].get_neighbors_within_tolerance(atom.crd(), 
						max_linker_length, tolerance);
				
				for (auto &patompoint : neighbors) {
					if (ptompoint->get_atom() == join_atom(segment1, segment2)) {
						clashed.insert(&patompoint->get_state());
					}
				}
			}
		}
		return join;
	}

	State::Set Poses::get_clashed_states(const State &state) const {
		State::Set clashed;
		for (auto &atom : state.get_atoms()) {
			AtomPoint::PVec neighbors = __grid.get_neighbors(atom.crd(), atom.radius());
			for (auto &patompoint : neighbors) {
				if (state.get_segment().get_id() != patompoint->get_state().get_segment().get_id()) {
					if (!excluded_atoms(atom, patompoint->get_atom(), state.get_segment(), patompoint->get_state().get_segment()))
						clashed.insert(&patompoint->get_state());
				}
			}
		}
		return clashed;
	}
};
