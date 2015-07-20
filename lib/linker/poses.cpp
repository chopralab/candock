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
	Poses::Poses(const Seed::Graph &seed_graph) {
		try {
			for (auto &seed : seed_graph) {
				auto &segment = seed.get_segment();
				for (auto &pstate : segment.get_states()) {
					auto &state = *pstate;
					for (int i = 0; i < state.get_crds().size(); ++i) {
						__atompoints[segment.get_id()].push_back(unique_ptr<AtomPoint>(new AtomPoint(state.get_crd(i), state.get_segment().get_atom(i), state)));
					}
				}
				__grid[segment.get_id()] = AtomPoint::Grid(__atompoints[segment.get_id()]);
			}
		} catch (...) {
			dbgmsg("FAILURE: constructor of Poses failed...");
			throw;
		}
	}
	
	Poses::Poses(const State::Set &states) {
		try {
			for (auto &pstate : states) {
				auto &state = *pstate;
				for (int i = 0; i < state.get_crds().size(); ++i) {
					__atompoints_small.push_back(unique_ptr<AtomPoint>(new AtomPoint(state.get_crd(i), state.get_segment().get_atom(i), state)));
				}
			}
			__grid_small = AtomPoint::Grid(__atompoints_small);
		} catch (...) {
			dbgmsg("FAILURE: constructor of Poses failed...");
			throw;
		}
	}
	
	State::Set Poses::get_clashed_states(const State &state, Segment &segment2) {
		State::Set clashed;
		AtomPoint::Grid &g = __grid.at(segment2.get_id());
		for (int i = 0; i < state.get_segment().get_atoms().size(); ++i) {
			AtomPoint::PVec neighbors = g.get_clashed(Molib::Atom(state.get_crd(i), state.get_segment().get_atom(i).idatm_type()));
			for (auto &patompoint : neighbors) {
				clashed.insert(&patompoint->get_state());
			}
		}
		return clashed;
	}

	State::Set Poses::get_clashed_states(const State &state) {
		State::Set clashed;
		AtomPoint::Grid &g = __grid_small;
		for (int i = 0; i < state.get_segment().get_atoms().size(); ++i) {
			AtomPoint::PVec neighbors = g.get_clashed(Molib::Atom(state.get_crd(i), state.get_segment().get_atom(i).idatm_type()));
			for (auto &patompoint : neighbors) {
				clashed.insert(&patompoint->get_state());
			}
		}
		return clashed;
	}

	State::Set Poses::get_join_states(const State &state, Segment &segment2, Molib::Atom::Pair &jatoms, const double max_linker_length, const double tol_seed_dist) {
		State::Set join;
		AtomPoint::Grid &g = __grid.at(segment2.get_id());
		for (int i = 0; i < state.get_segment().get_atoms().size(); ++i) {
			if (&state.get_segment().get_atom(i) == jatoms.first) {
				AtomPoint::PVec neighbors = g.get_neighbors_within_tolerance(state.get_crd(i), 
					max_linker_length, tol_seed_dist);
						
				for (auto &patompoint : neighbors) {
					if (&patompoint->get_atom() == jatoms.second) {
						join.insert(&patompoint->get_state());
					}
				}
			}
		}
		return join;
	}
};
