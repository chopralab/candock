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
	Poses::Poses(const Seed::Graph &seed_graph, const double tol_seed_dist) : __tol_seed_dist(tol_seed_dist) {
		try {
			__initialize_grid(seed_graph);
		} catch (...) {
			dbgmsg("FAILURE: constructor of Poses failed...");
			throw;
		}
	}
	
	void Poses::__initialize_grid(const Seed::Graph &seed_graph) {
		for (auto &seed : seed_graph) {
			auto &segment = seed.get_segment();
			vector<AtomPoint*> ap;
			for (auto &pstate : segment.get_states()) {
				auto &state = *pstate;
				for (auto &kv : state.get_atoms()) {
					Molib::Atom &atom = const_cast<Molib::Atom&>(*kv.first);
					__atompoints.push_back(unique_ptr<AtomPoint>(new AtomPoint(atom, state)));
					ap.push_back(&*__atompoints.back());
				}
			}
			__grid[segment.get_id()] = Grid<AtomPoint>(ap);
		}
	}
	
	State::Set Poses::get_clashed_states(const State &state, Segment &segment2) {
		State::Set clashed;
		Grid<AtomPoint> &g = __grid.at(segment2.get_id());
		for (auto &kv : state.get_atoms()) {
			Molib::Atom &atom = const_cast<Molib::Atom&>(*kv.first);
			AtomPoint::PVec neighbors = g.get_neighbors(atom.crd(), atom.radius());
			for (auto &patompoint : neighbors) {
				clashed.insert(&patompoint->get_state());
			}
		}
		return clashed;
	}

	State::Set Poses::get_join_states(const State &state, Segment &segment2, Molib::Atom::Pair &jatoms, const double max_linker_length) {
		State::Set join;
		Grid<AtomPoint> &g = __grid.at(segment2.get_id());
		for (auto &kv : state.get_atoms()) {
			Molib::Atom &atom = const_cast<Molib::Atom&>(*kv.first);
			if (&atom == jatoms.first) {
				AtomPoint::PVec neighbors = g.get_neighbors_within_tolerance(atom.crd(), 
					max_linker_length, __tol_seed_dist);
						
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
