/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "poses.hpp"
#include "state.hpp"
#include "segment.hpp"
#include "molib/molecule.hpp"
#include "molib/bond.hpp"
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
