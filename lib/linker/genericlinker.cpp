#include "linker.hpp"
#include "score/score.hpp"
#include "molib/nrset.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "modeler/modeler.hpp"
#include "geom3d/quaternion.hpp"
#include "poses.hpp"
#include "molib/internal.hpp"
#include <queue>

using namespace std;

namespace Linker {

	/**
	 * Each state is a mapping of docked seed atoms to ligands's 
	 * segment atoms
	 */
	void Linker::GenericLinker::__create_states(const Segment::Graph &segment_graph, const Molib::NRset &top_seeds) {
		for (auto &seed_mols : top_seeds) {
			const int seed_id = stoi(seed_mols.name());
			dbgmsg("seed id = " << seed_id);
			// loop over segments with this name, seeds coordinates will not change
			for (auto &segment : segment_graph) {
				if (segment.get_seed_id() == seed_id) {

					dbgmsg("segment seed id = " << segment.get_seed_id());
					// create a graph out of the seed "nm" of "to be" ligand molecule
					Molib::Atom::Graph gs = Molib::Atom::create_graph(seed_mols.first().get_atoms());
					dbgmsg("gs = " << endl << gs);
					Molib::Atom::Graph gd = Molib::Atom::create_graph(segment.get_atoms());
					dbgmsg("gd = " << endl << gd);
					// map seed molecules atom numbers to ligand molecule
					Molib::Atom::Graph::Matches m = gd.match(gs);
					dbgmsg("__create_states : m.size() = " << m.size());
					// consider ONLY THE FIRST mapping of seed to seed (the symmetry has 
					// been taken care of in the rigid docking itself)...
					if (m.size() < 1) throw Error ("die : could not find mapping of docked seed to segment with same name ???");
					auto &mv = *m.begin();
					// create a state
					auto &vertices1 = mv.first;
					auto &vertices2 = mv.second;
					// for every docked rigid fragment
					for (auto &seed_molecule : seed_mols) {
						Geom3D::Point::Vec crds(vertices2.size());
						const Molib::Atom::Vec &seed_atoms = seed_molecule.get_atoms();
						for (size_t i = 0; i < vertices2.size(); ++i) {
							crds[vertices1[i]] = seed_atoms.at(vertices2[i])->crd();
							dbgmsg("adding matched vertex pair " << vertices1[i] 
								<< "," << vertices2[i] << " new coordinates of atom " 
								<< segment.get_atom(vertices1[i]).atom_number() << " are " 
								<< crds[vertices1[i]]);
						}
						const double energy = stod(seed_molecule.name());
						if (energy < __max_allow_energy) {
							dbgmsg("adding docked state " << segment.get_seed_id() << " with energy of " << energy);
							// ONLY COPY SEGMENT COORDS NOT SEED
							segment.add_state(unique_ptr<State>(new State(segment, crds, energy)));
						}
					}
				}
			}
		}
	}
	
	/**
	 * Initialize OpenMM: add ligand coordinates (random) since otherwise 
	 * there are zeroes (at iterative minimization, when masking out parts 
	 * of ligand) in the openmm positions array which openmm does not like
	 * 
	 */
	void Linker::GenericLinker::init_openmm() {
		__modeler.add_topology(__receptor.get_atoms());
		__modeler.add_topology(__ligand.get_atoms());
		
		__modeler.init_openmm();

		__modeler.add_random_crds(__ligand.get_atoms());
	}
	
	Partial::Vec Linker::GenericLinker::init_conformations() {
		log_step << "Starting connection of seeds for ligand " << __ligand.name() << endl;
		
		// A graph of segments is constructed in which each rigid segment is 
		// a vertex & segments that are part of seeds are identified.		
		__segment_graph = Segment::create_graph(__ligand);
		if (!__segment_graph.find_cycles_connected_graph().empty()) {
			throw Error("die : cyclic molecules are currently not supported");
		}
		
		dbgmsg("segment graph for ligand " << __ligand.name() << " = " << __segment_graph);
		__create_states(__segment_graph, __top_seeds);
		
		__seed_graph = Seed::create_graph(__segment_graph);
		dbgmsg("seed graph for ligand " << __ligand.name() << " = " << __seed_graph);
		
		return __generate_rigid_conformations(__seed_graph);
	}
	
	/**
	 * Link the conformations using the A-STAR algorithm
	 * 
	 */
	DockedConformation::Vec Linker::GenericLinker::compute_conformations(const Partial::Vec &partials) {

		DockedConformation::Vec docked_conformations;

		size_t failed_connections = 0;

		for (auto &partial : partials) {
			auto &conformation = partial.get_states();
				
			if (!__has_blacklisted(conformation, __blacklist)) {
				try {
					Partial grow_link_energy(conformation, partial.get_energy());
					dbgmsg("connecting ligand " << __ligand.name() 
						<< endl << "possible starting conformation is " << endl
						<< grow_link_energy);
					vector<unique_ptr<State>> states;
					docked_conformations.push_back(
						__a_star(__segment_graph.size(), grow_link_energy, states, __link_iter));
				} catch (ConnectionError& e) {
					dbgmsg("exception : " << e.what());
					++failed_connections;
					for (auto &failed_pair : e.get_failed_state_pairs())
						__blacklist.insert(failed_pair);
				}
			}
		}

		if (failed_connections >= partials.size())
			throw Error ("No generated confirmations are valid.");

		if (__max_iterations_final > 0)
			return __minimize_final(docked_conformations);

                // Potenital energy will never be set, and that's ok
		for (auto &docked : docked_conformations)
			docked.get_receptor().undo_mm_specific();
		
		return docked_conformations;
	}

	/**
	 * Final minimization of each docked ligand conformation
	 * with full ligand and receptor flexibility
	 * 
	 */
	DockedConformation::Vec Linker::GenericLinker::__minimize_final(DockedConformation::Vec &docked_conformations) {

		DockedConformation::Vec minimized_conformations;

		for (auto &docked : docked_conformations) {
		
			try {

				__modeler.add_crds(__receptor.get_atoms(), docked.get_receptor().get_crds());
				__modeler.add_crds(__ligand.get_atoms(), docked.get_ligand().get_crds());
				
				__modeler.init_openmm_positions();
				
				__modeler.unmask(__receptor.get_atoms());
				__modeler.unmask(__ligand.get_atoms());
		
				__modeler.set_max_iterations(__max_iterations_final); // until converged
				__modeler.minimize_state(__ligand, __receptor, __score);

                                __modeler.mask(__receptor.get_atoms());
                                const double potential_energy = __modeler.potential_energy();

				// init with minimized coordinates
				Molib::Molecule minimized_receptor(__receptor, __modeler.get_state(__receptor.get_atoms()));
				Molib::Molecule minimized_ligand(__ligand, __modeler.get_state(__ligand.get_atoms()));
		
				minimized_receptor.undo_mm_specific();
				
				Molib::Atom::Grid gridrec(minimized_receptor.get_atoms());
		
				minimized_conformations.push_back(DockedConformation(minimized_ligand, minimized_receptor,
					__score.non_bonded_energy(gridrec, minimized_ligand), potential_energy));
		
			} catch(OMMIface::Modeler::MinimizationError &e) {
				log_error << "MinimizationError: skipping minimization of one conformation of ligand " 
					<< __ligand.name() << " due to : " << e.what() << endl;
			}
		}
		return minimized_conformations;
	}
	
	bool Linker::GenericLinker::__clashes_receptor(const State &current) const {
		for (size_t i = 0; i < current.get_crds().size(); ++i) {
			const Molib::Atom &a = current.get_segment().get_atom(i); 
			const Geom3D::Coordinate &c = current.get_crd(i); 
			dbgmsg("in clashes_receptor test coordinate = " << c);
			if (__gridrec.clashes(Molib::Atom(c, a.idatm_type()), __clash_coeff)) return true;
		}
		return false;
	}
	
	bool Linker::GenericLinker::__clashes_ligand(const State &current, 
		const Partial &conformation, const State &) const {

		// clashes between current segment and previous segments
		for (auto &pstate : conformation.get_states()) {
			dbgmsg("clashes_ligand test for state " << current 
				<< " and state " << *pstate);
			if (current.clashes(*pstate, __clash_coeff)) 
				return true;
		}
		return false;
	}

	double Linker::GenericLinker::__distance(const State &start, const State &goal) const {
		const Segment &start_segment = start.get_segment();
		const Segment &goal_segment = goal.get_segment();
		const int start_atom_idx = start_segment
			.get_bond(start_segment.get_next(goal_segment)).idx1();
		const int goal_atom_idx = goal_segment
			.get_bond(goal_segment.get_next(start_segment)).idx1();
		const Geom3D::Coordinate &first_crd = start.get_crd(start_atom_idx);
		const Geom3D::Coordinate &last_crd = goal.get_crd(goal_atom_idx);
		return first_crd.distance(last_crd);
	}
	
	Geom3D::Point::Vec Linker::GenericLinker::__rotate(const Geom3D::Quaternion &q, 
		const Geom3D::Point &p1, const Geom3D::Point::Vec &crds) {

		Geom3D::Point::Vec rotated;
		for (auto &crd : crds) {	
			rotated.push_back(q.rotatedVector(crd - p1) + p1); 
		}
		return rotated;
	}
	
	State::Vec Linker::GenericLinker::__compute_neighbors(const State &curr_state, Segment &next,
		vector<unique_ptr<State>> &states) {
		const int ini_sz = states.size();
		dbgmsg("in compute_neighbors for state = " << curr_state);
		const Segment &current = curr_state.get_segment();
		dbgmsg("compute_neighbors : current segment is " << current);
		dbgmsg("compute_neighbors : next segment is " << next);
		const Molib::Bond &btorsion = current.get_bond(next);
		const int idx3 = btorsion.idx2();
		const Molib::Atom *a3 = &btorsion.atom2();
		const int idx2 = btorsion.idx1(); 
		const Molib::Atom *a2 = &btorsion.atom1();
		const int idx1 = current.adjacent_in_segment(*a2, *a3);	// segments overlap on 1 rotatable bond
		const Molib::Atom *a1 = &current.get_atom(idx1);
		const int idx4 = next.adjacent_in_segment(*a3, *a2);
		const Molib::Atom *a4 = &next.get_atom(idx4);
		
		dbgmsg("a4 = " << *a4 << " coordinates not set yet!");
		Geom3D::Coordinate crd3(curr_state.get_crd(idx3));
		dbgmsg("a3 = " << *a3 << " crd3 = " << curr_state.get_crd(idx3));
		Geom3D::Coordinate crd2(curr_state.get_crd(idx2));
		dbgmsg("a2 = " << *a2 << " crd2 = " << curr_state.get_crd(idx2));
		Geom3D::Coordinate crd1(curr_state.get_crd(idx1));
		dbgmsg("a1 = " << *a1 << " crd1 = " << curr_state.get_crd(idx1));
		
		const double saved_dihedral = __ic.get_dihedral(*a1, *a2, *a3, *a4);
		__ic.set_dihedral(*a1, *a2, *a3, *a4, 0.0);
		
		states.push_back(unique_ptr<State>(new State(next, 
			__ic.cartesian(*a1, *a2, *a3, crd1, crd2, crd3, next.get_atoms()))));

		__ic.set_dihedral(*a1, *a2, *a3, *a4, saved_dihedral);

#ifndef NDEBUG
		State &initial = *states.back();
#endif
		dbgmsg("this is initial state = " << initial);
		dbgmsg("a4 = " << *a4 << " newly determined crd4 = " 
			<< initial.get_crd(idx4));
		dbgmsg("rotate next segment on vector = "
				<< crd2 << " - " << crd3
				<< " by " << Geom3D::degrees(__spin_degrees) 
				<< " degree increments");
		const Geom3D::Quaternion q(Geom3D::Vector3(crd3 - crd2).norm()*sin(__spin_degrees), cos(__spin_degrees));
		for (double angle = __spin_degrees; angle < M_PI; angle += __spin_degrees) {
			State &previous_rotated = *states.back();
			states.push_back(unique_ptr<State>(new State(next, 
				__rotate(q, crd2, previous_rotated.get_crds()))));
			dbgmsg("rotated state at angle = " << Geom3D::degrees(angle)
				<< " is " << *states.back());
		}
		
		State::Vec ret;
		for (size_t i = ini_sz; i < states.size(); ++i) ret.push_back(&*states[i]);
	
		return ret;
	}
	
	bool Linker::GenericLinker::__check_distances_to_seeds(const State &curr_state, 
		const Segment &adjacent, const SegStateMap &docked_seeds) {

		const Segment &current = curr_state.get_segment();
		dbgmsg("curr_state= " << curr_state);
		dbgmsg("adjacent= " << adjacent);
		for (auto &adjacent_seed_segment : adjacent.get_adjacent_seed_segments()) {
			auto it = docked_seeds.find(adjacent_seed_segment);
			if (it != docked_seeds.end() && it->first != &current) {

				const double distance = __distance(curr_state, *it->second);
				const double mll = current.get_max_linker_length(*it->first);

				dbgmsg("distance = " << distance);
				dbgmsg("mll = " << mll);
				dbgmsg("is distance OK ? " << boolalpha << (distance < mll + __tol_seed_dist
					&& distance > mll - __tol_seed_dist));
				dbgmsg("curr_state= " << curr_state);
				dbgmsg("adjacent seed state= " << *it->second);

				if (distance > mll + __tol_seed_dist || distance < mll - __tol_seed_dist) {
					dbgmsg("returning false");
					return false;
				}
			}
		}
		dbgmsg("returning true");
		return true;
	}

	State* Linker::GenericLinker::__is_seed(const Segment &seg, const SegStateMap &docked_seeds) {
		auto it = docked_seeds.find(&seg);
		return (it == docked_seeds.end() ? nullptr : it->second);
	}

	string Linker::GenericLinker::to_pdb(const Partial &conformation) {
		stringstream os;
		os << "MODEL" << endl;
		os << "REMARK   8 ENERGY " << conformation.get_energy() << endl;
		for (auto &pstate : conformation.get_states())
			os << pstate->pdb();
		os << "ENDMDL" << endl;
		return os.str();
	}

	pair<State*, Segment*> Linker::GenericLinker::__find_good_neighbor(const Partial &curr_conformation, const SegStateMap &docked_seeds) {

		Segment::Set done_segments, free_seeds;
		for (auto &pcurr_state : curr_conformation.get_states()) {
			done_segments.insert(const_cast<Segment*>(&pcurr_state->get_segment()));
		}
		for (auto &kv : docked_seeds) {
			const Segment &docked_segment = *kv.first;
			if (!done_segments.count(const_cast<Segment*>(&docked_segment))) {
				free_seeds.insert(const_cast<Segment*>(&docked_segment));
			}
		}
		pair<State*, Segment*> good;
		for (auto &pcurr_state : curr_conformation.get_states()) {
			State &curr_state = *pcurr_state;
			const Segment &curr_segment = curr_state.get_segment();
			for (auto &pgoal : free_seeds) {
				if (curr_segment.has_next(*pgoal)) {
					Segment &next = curr_segment.get_next(*pgoal);
					if (!done_segments.count(&next)) {
						return {&curr_state, &next};
					}
				}
			}
		}
		for (auto &pcurr_state : curr_conformation.get_states()) {
			State &curr_state = *pcurr_state;
			for (auto &adj : curr_state.get_segment()) { 
				if (!done_segments.count(&adj)) { // don't go back
					return {&curr_state, &adj};
				}
			}
		}

		throw Error ("Could not find a good neighbor");
	}

	Array2d<bool> Linker::GenericLinker::__find_compatible_state_pairs(const Seed::Graph &seed_graph, const int sz) {
		/* Find all pairs of compatible states at correct distances for 
		 * the multi-seed molecules
		 */
		
		Array2d<bool> conn(sz);
		Poses poses(seed_graph);
		for (size_t u = 0; u < seed_graph.size(); ++u) {
			Segment &segment1 = seed_graph[u].get_segment();
			for (size_t v = u + 1; v < seed_graph.size(); ++v) {
				Segment &segment2 = seed_graph[v].get_segment();
				const double max_linker_length = segment1.get_max_linker_length(segment2);
				Molib::Atom::Pair jatoms{&segment1.get_bond(segment1.get_next(segment2)).atom1(), 
					&segment2.get_bond(segment2.get_next(segment1)).atom1()};
#ifndef NDEBUG
				const Molib::Bond &excluded = (segment1.is_adjacent(segment2) 
					? segment1.get_bond(segment2) : Molib::Bond());
				dbgmsg("finding compatible state pairs for segments " << segment1 << " and " 
					<< segment2 << " with maximum linker length " << max_linker_length << " and "
					<< " join atom1 = " << jatoms.first->atom_number() << " join atom2 = "
					<< jatoms.second->atom_number() << " excluded = " << excluded);
#endif
				for (auto &pstate1 : segment1.get_states()) {
					State &state1 = *pstate1;
					State::Set join_states = poses.get_join_states(state1, segment2, jatoms, 
						max_linker_length, __lower_tol_seed_dist, __upper_tol_seed_dist);

					const int i = state1.get_id();
					dbgmsg("comparing state1(" << i << ") = " << state1);
					for (auto &pstate2 : join_states) {
						dbgmsg("with state2(" << pstate2->get_id() << ") = " << *pstate2);
						if (!state1.clashes(*pstate2, __clash_coeff)) {
							const int j = pstate2->get_id();
							conn.set(i, j);
							conn.set(j, i);
							dbgmsg("compatible states " << i << " and " << j);
						}
					}
				}
			}
		}

		return conn;
	}

	bool Linker::GenericLinker::__has_blacklisted(const State::Vec &conformation, 
		const set<State::ConstPair> &blacklist) {

		for (size_t i = 0; i < conformation.size(); ++i) {
			for (size_t j = i + 1; j < conformation.size(); ++j) {
				if (blacklist.count({conformation[i], conformation[j]})
					|| blacklist.count({conformation[j], conformation[i]}))
					return true;
			}
		}
		return false;
	}

};
