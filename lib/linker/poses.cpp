#include "candock/linker/poses.hpp"
#include "candock/linker/state.hpp"
#include "candock/linker/segment.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/molib/bond.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/help.hpp"
#include "candock/helper/array2d.hpp"
#include "candock/graph/mcqd.hpp"
#include <queue>

using namespace std;

namespace Linker {
	Poses::Poses(const Seed::Graph &seed_graph) {
		try {
			for (auto &seed : seed_graph) {
				auto &segment = seed.get_segment();
				for (auto &pstate : segment.get_states()) {
					auto &state = *pstate;
					for (size_t i = 0; i < state.get_crds().size(); ++i) {
						if (segment.is_join_atom(i)) //FIXME What is meant here? Technically only the debug statement is run, but is this check supposed to run __atompoints too?
							dbgmsg("pushing to atompoints crd = " << state.get_crd(i) << " atom = " << segment.get_atom(i) << " for STATE = " << state);
							__atompoints[segment.get_id()].push_back(unique_ptr<AtomPoint>(new AtomPoint(state.get_crd(i), segment.get_atom(i), state)));
					}
				}
				__grid[segment.get_id()] = AtomPoint::Grid(__atompoints[segment.get_id()]);
			}
		} catch (...) {
			dbgmsg("FAILURE: constructor of Poses failed...");
			throw;
		}
	}
	
	State::Set Poses::get_join_states(const State &state, Segment &segment2, Molib::Atom::Pair &jatoms, 
		const double max_linker_length, const double lower_tol_seed_dist, const double upper_tol_seed_dist) {
		State::Set join;
		AtomPoint::Grid &g = __grid.at(segment2.get_id());
		for (size_t i = 0; i < state.get_segment().get_atoms().size(); ++i) {
			dbgmsg("atom1 = " << state.get_segment().get_atom(i).atom_number() << " atom2 = " << jatoms.first->atom_number()
				<< " equal = " << boolalpha << (&state.get_segment().get_atom(i) == jatoms.first));
			if (&state.get_segment().get_atom(i) == jatoms.first) {
				dbgmsg("crd = " << state.get_crd(i) << " mll = " << max_linker_length 
					<< " lower_tol_seed_dist = " << lower_tol_seed_dist 
					<< " upper_tol_seed_dist = " << upper_tol_seed_dist); 
				AtomPoint::PVec neighbors = g.get_neighbors_within_tolerance_asymmetric(state.get_crd(i), 
					max_linker_length, lower_tol_seed_dist, upper_tol_seed_dist);
				dbgmsg("found neighbors = " << neighbors.size());
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
