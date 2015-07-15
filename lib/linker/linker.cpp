#include "linker.hpp"
#include "geom3d/quaternion.hpp"
#include "score/score.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "graph/mcqd.hpp"
#include <queue>

using namespace std;

namespace Molib {
	ostream& operator<<(ostream& os, const Linker::Conf &conf)	{
		os << "start link ++++++++++++++++++++++++++++++" << endl;
		os << "ENERGY = " << conf.second << endl;
		for (auto &state : conf.first)
			os << state << endl;
		os << "end link --------------------------------" << endl;
		return os;
	}

	ostream& operator<<(ostream& os, const Linker::LinkEnergy &le)	{
		os << "start link ++++++++++++++++++++++++++++++" << endl;
		os << "ENERGY = " << le.second << endl;
		for (auto &pstate : le.first)
			os << *pstate << endl;
		os << "end link --------------------------------" << endl;
		return os;
	}

	ostream& operator<<(ostream& os, const vector<Linker::LinkEnergy> &vec_le)	{
		for (auto &le : vec_le) {
			os << le << endl;
		}			
		return os;
	}

	bool Linker::__link_adjacent(const Segment::Graph::Path &path) {
		/* Returns true if path has no seed segments along the way
		 * 
		 */
#ifndef NDEBUG
		for (auto it = path.begin(); it != path.end(); ++it) {
			dbgmsg(**it << " is seed = " << boolalpha << (*it)->is_seed());
		}
#endif
		for (auto it = path.begin() + 1; it != path.end() - 1; ++it) {
			dbgmsg(**it);
			if ((*it)->is_seed())
				return false;
		}
		return true;
	}

	void Linker::__init_max_linker_length(const Segment::Paths &paths) {
		for (auto &kv : paths) {
			auto &seg_pair = kv.first;
			Segment::Graph::Path path(kv.second.begin(), kv.second.end());
			__compute_max_linker_length(path);
			reverse(path.begin(), path.end());
			__compute_max_linker_length(path);
		}
	}

	void Linker::__compute_max_linker_length(Segment::Graph::Path &path) {
		double d = 0.0;
		int i = 0;
		for (; i < path.size() - 2; i += 2) {
			const Bond &b1 = path[i]->get_bond(*path[i + 1]);
			const Bond &b2 = path[i + 1]->get_bond(*path[i + 2]);
			path[0]->set_max_linker_length(*path[i + 1], d + b1.length());
			path[i + 1]->set_max_linker_length(*path[0], d + b1.length());
			//~ d += b1.first_atom().crd().distance(b2.second_atom().crd());
			d += b1.atom1().crd().distance(b2.atom2().crd());
			path[0]->set_max_linker_length(*path[i + 2], d);
			path[i + 2]->set_max_linker_length(*path[0], d);
			dbgmsg("max_linker_length between " << *path[0] << " and " 
				<< *path[i + 1] << " = " << path[0]->get_max_linker_length(*path[i + 1]));
			dbgmsg("max_linker_length between " << *path[0] << " and " 
				<< *path[i + 2] << " = " << path[0]->get_max_linker_length(*path[i + 2]));
		}
		if (i < path.size() - 1) {
			const Bond &b = path[i]->get_bond(*path[i + 1]);
			d += b.length();
			path[0]->set_max_linker_length(*path[i + 1], d);
			path[i + 1]->set_max_linker_length(*path[0], d);
			dbgmsg("last max_linker_length between " << *path[0] << " and " 
				<< *path[i + 1] << " = " << path[0]->get_max_linker_length(*path[i + 1]));
		}
		dbgmsg("TOTAL max_linker_length between " << *path[0] << " and " 
			<< *path[path.size() - 1] << " = " 
			<< path[0]->get_max_linker_length(*path[path.size() - 1]));
	}

	void Linker::__create_states(const Segment::Graph &segment_graph, const NRset &top_seeds) {
		/* each state is a mapping of docked seed atoms to ligands's 
		 * segment atoms
		 */
		for (auto &seed_mols : top_seeds) {
			const int seed_id = stoi(seed_mols.name());
			dbgmsg("seed id = " << seed_id);
			// loop over segments with this name, seeds coordinates will not change
			for (auto &segment : segment_graph) {
				if (segment.get_seed_id() == seed_id) {

					dbgmsg("segment seed id = " << segment.get_seed_id());
					// create a graph out of the seed "nm" of "to be" ligand molecule
					MolGraph gs = create_graph(seed_mols.first().get_atoms());
					dbgmsg("gs = " << endl << gs);
					MolGraph gd = create_graph(segment.get_atoms());
					dbgmsg("gd = " << endl << gd);
					// map seed molecules atom numbers to ligand molecule
					MolGraph::Matches m = gd.match(gs);
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
						//~ AtomVec ordered_atoms;
						AtomToCrd atom_crd;
						for (int i = 0; i < vertices2.size(); ++i) {
							Atom &v1 = gd[vertices1[i]];
							Atom &v2 = *seed_molecule.get_atoms()[vertices2[i]];
							atom_crd[&v1] = v2.crd();
							
							dbgmsg("adding matched vertex pair " << vertices1[i] 
								<< "," << vertices2[i] << " new coordinates of atom " 
								<< v1.atom_number() << " are " 
								<< v2.crd());
						}
						const double energy = stod(seed_molecule.name());
						dbgmsg("adding docked state " << segment.get_seed_id() << " with energy of " << energy);
						// ONLY COPY SEGMENT COORDS NOT SEED
						segment.add_state(unique_ptr<State>(new State(segment, atom_crd, energy)));
					}
				}
			}
		}
	}
	
	static double torsion_energy(const State &first, const State &second) { return 0.0; }

	bool Linker::__clashes_receptor(const State &current) const {
		for (auto &kv : current.get_atoms()) {
			const Atom &a = *kv.first; 
			const Geom3D::Coordinate &c = kv.second; 
			dbgmsg("in clashes_receptor test coordinate = " << c);
			if (__gridrec.clashes(Atom(c, a.idatm_type()))) return true;
		}
		return false;
	}
	bool Linker::__clashes_ligand(const State &current, 
		const LinkEnergy &conformation, const State &prev) const {

		const Bond &excluded = current.get_segment().get_bond(prev.get_segment());
		// clashes between current segment and previous segments
		for (auto &pstate : conformation.first) {
			dbgmsg("clashes_ligand test for state " << current 
				<< " and state " << *pstate);
			if (current.clashes(*pstate, excluded)) 
				return true;
		}
		return false;
	}

	double Linker::__distance(const State &start, const State &goal) const {
		const Segment &start_segment = start.get_segment();
		const Segment &goal_segment = goal.get_segment();
		const Atom &first_atom = start_segment
			//~ .get_bond(start_segment.get_next(goal_segment)).first_atom();
			.get_bond(start_segment.get_next(goal_segment)).atom1();
		const Atom &last_atom = goal_segment
			//~ .get_bond(goal_segment.get_next(start_segment)).first_atom();
			.get_bond(goal_segment.get_next(start_segment)).atom1();
		const Geom3D::Coordinate &first_crd = start.get_atom_crd(first_atom);
		const Geom3D::Coordinate &last_crd = goal.get_atom_crd(last_atom);
		return first_crd.distance(last_crd);
	}
	
	AtomToCrd Linker::__rotate(const Geom3D::Quaternion &q, 
		const Geom3D::Point &p1, const Geom3D::Point &p2, const AtomToCrd &atom_crd) {
		AtomToCrd new_crd;
		for (auto &kv : atom_crd) {	
			new_crd.insert({kv.first, q.rotatedVector(kv.second - p1) + p1}); 
		}
		return new_crd;
	}
	
	State::Vec Linker::__compute_neighbors(const State &curr_state, Segment &next,
		vector<unique_ptr<State>> &states) {
		dbgmsg("in compute_neighbors for state = " << curr_state);
		State::Vec l;
		const Segment &current = curr_state.get_segment();
		dbgmsg("compute_neighbors : current segment is " << current);
		dbgmsg("compute_neighbors : next segment is " << next);
		const Bond &btorsion = current.get_bond(next);
		//~ const Atom *a3 = &btorsion.second_atom(); 
		const Atom *a3 = &btorsion.atom2(); 
		//~ const Atom *a2 = &btorsion.first_atom();
		const Atom *a2 = &btorsion.atom1();
		const Atom *a1 = &current.adjacent_in_segment(*a2, *a3);	// segments overlap on 1 rotatable bond
		const Atom *a4 = &next.adjacent_in_segment(*a3, *a2);
		dbgmsg("a4 = " << *a4 << " coordinates not set yet!");
		Geom3D::Coordinate crd3(curr_state.get_atom_crd(*a3));
		dbgmsg("a3 = " << *a3 << " crd3 = " << curr_state.get_atom_crd(*a3));
		Geom3D::Coordinate crd2(curr_state.get_atom_crd(*a2));
		dbgmsg("a2 = " << *a2 << " crd2 = " << curr_state.get_atom_crd(*a2));
		Geom3D::Coordinate crd1(curr_state.get_atom_crd(*a1));
		dbgmsg("a1 = " << *a1 << " crd1 = " << curr_state.get_atom_crd(*a1));
		states.push_back(unique_ptr<State>(new State(next, 
			__ic.cartesian(*a1, *a2, *a3, crd1, crd2, crd3, next.get_atoms()))));
		State &initial = *states.back();
		dbgmsg("this is initial state = " << initial);
		dbgmsg("a4 = " << *a4 << " newly determined crd4 = " 
			<< Geom3D::Coordinate(initial.get_atom_crd(*a4)));
		l.push_back(&initial);
		dbgmsg("rotate next segment on vector = "
				<< crd2 << " - " << crd3
				<< " by " << Geom3D::degrees(__spin_degrees) 
				<< " degree increments");
		const Geom3D::Quaternion q(Geom3D::Vector3(crd3 - crd2).norm()*sin(__spin_degrees), cos(__spin_degrees));
		for (double angle = __spin_degrees; angle < M_PI; angle += __spin_degrees) {
			State &previous_rotated = *states.back();
			states.push_back(unique_ptr<State>(new State(next, 
				__rotate(q, crd2, crd3, previous_rotated.get_atoms()))));
			dbgmsg("rotated state at angle = " << Geom3D::degrees(angle)
				<< " is " << *states.back());
			l.push_back(&*states.back());
		}
		return l;
	}
	
	bool Linker::__check_distances_to_seeds(const State &curr_state, 
		const Segment &adjacent, const SegStateMap &docked_seeds) {
		const double mll_eps = 0.001; // add to max linker length to avoid
		// checks to "back" seeds to be true, e.g., 2.44893690 > 2.44893690 
		// should be false (but sometimes it returns true)
		const Segment &current = curr_state.get_segment();
		dbgmsg("curr_state= " << curr_state);
		dbgmsg("adjacent= " << adjacent);
		for (auto &adjacent_seed_segment : adjacent.get_adjacent_seed_segments()) {
			auto it = docked_seeds.find(adjacent_seed_segment);
			if (it != docked_seeds.end() && it->first != &current) {
				dbgmsg("distance = " << __distance(curr_state, *it->second));
				dbgmsg("mll = " << current.get_max_linker_length(*it->first));
				dbgmsg("distance > mll = " << boolalpha << (__distance(curr_state, *it->second) 
					> current.get_max_linker_length(*it->first) + mll_eps));
				dbgmsg("curr_state= " << curr_state);
				dbgmsg("adjacent seed state= " << *it->second);
				if (__distance(curr_state, *it->second) > current
					.get_max_linker_length(*it->first) + mll_eps) {
					dbgmsg("returning false");
					return false;
				}
			}
		}
		dbgmsg("returning true");
		return true;
	}

	State* Linker::__is_seed(const Segment &seg, const SegStateMap &docked_seeds) {
		auto it = docked_seeds.find(&seg);
		return (it == docked_seeds.end() ? nullptr : it->second);
	}

	string Linker::to_pdb(const LinkEnergy &conformation) {
		stringstream os;
		os << "MODEL" << endl;
		os << "REMARK   8 ENERGY " << conformation.second << endl;
		for (auto &pstate : conformation.first)
			os << pstate->pdb();
		os << "ENDMDL" << endl;
		return os.str();
	}

	pair<State*, Segment*> Linker::__find_good_neighbor(
		const LinkEnergy &curr_conformation, const SegStateMap &docked_seeds) {

		Segment::Set done_segments, free_seeds;
		for (auto &pcurr_state : curr_conformation.first) {
			done_segments.insert(const_cast<Segment*>(&pcurr_state->get_segment()));
		}
		for (auto &kv : docked_seeds) {
			const Segment &docked_segment = *kv.first;
			if (!done_segments.count(const_cast<Segment*>(&docked_segment))) {
				free_seeds.insert(const_cast<Segment*>(&docked_segment));
			}
		}
		pair<State*, Segment*> good;
		for (auto &pcurr_state : curr_conformation.first) {
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
		for (auto &pcurr_state : curr_conformation.first) {
			State &curr_state = *pcurr_state;
			for (auto &adj : curr_state.get_segment()) { 
				if (!done_segments.count(&adj)) { // don't go back
					return {&curr_state, &adj};
				}
			}
		}
	}

	Linker::Conf Linker::__a_star(const int segment_graph_size, 
		const LinkEnergy &start_conformation, int iter) {

		vector<unique_ptr<State>> states;
		for (auto &pstate : start_conformation.first)
			states.push_back(unique_ptr<State>(new State(*pstate)));
		if (start_conformation.first.empty())
			throw Error ("die : at least one docked anchor state is required for linking");
		SegStateMap docked_seeds;
		for (auto &pstate : states) 
			docked_seeds.insert({&pstate->get_segment(), &*pstate});
		set<State::ConstPair>	failed_state_pairs;
		LinkEnergy min_conformation;
		double min_energy = MAX_ENERGY;
		PriorityQueue openset; // openset has conformations sorted from lowest to highest energy
		openset.insert(LinkEnergy{State::Vec{&*states[0]}, states[0]->get_energy()});
		while(!openset.empty()) {
			if (--iter < 0) break;
			LinkEnergy curr_conformation = *openset.begin();
			openset.erase(openset.begin());
			dbgmsg("openset.size() = " << openset.size());
			dbgmsg("curr_conformation at step = " 
				<< iter << " = " << endl << curr_conformation
				<< endl << to_pdb(curr_conformation));
			if (curr_conformation.first.size() == segment_graph_size) {
				dbgmsg("CANDIDATE for minimum energy conformation at step = " 
					<< iter << " = " << endl << curr_conformation);
				const double curr_energy = curr_conformation.second;
				if (curr_energy < min_energy) {
					min_conformation = curr_conformation;
					min_energy = curr_energy;
					dbgmsg("ACCEPTED");
				}
			} else {
				// grow in the direction of (a) first "free" seed state (b) if such
				// path does not exist, then grow in the first possible direction
				pair<State*, Segment*> ret = __find_good_neighbor(curr_conformation,
					docked_seeds);
				State &curr_state = *ret.first;
				Segment &adj = *ret.second;
				State* adj_is_seed = __is_seed(adj, docked_seeds);
				// check seed distances here 
				dbgmsg("adj_is_seed " << (adj_is_seed ? "true" : "false"));
				dbgmsg("check_distances_to_seeds = " << boolalpha 
					<< (adj_is_seed || __check_distances_to_seeds(curr_state, adj, docked_seeds)));
				if (adj_is_seed || __check_distances_to_seeds(curr_state, adj, docked_seeds)) {
					auto ret2 = (adj_is_seed ? State::Vec{adj_is_seed} : 
						__compute_neighbors(curr_state, adj, states));
					for (auto &pneighbor : ret2) { 
						dbgmsg("CHECKING NEIGHBOR : " << *pneighbor);
						dbgmsg("clashes_receptor = " << boolalpha 
							<< __clashes_receptor(*pneighbor));
						dbgmsg("clashes_ligand = " << boolalpha 
							<< __clashes_ligand(*pneighbor, curr_conformation, curr_state));
						if (!__clashes_receptor(*pneighbor) 
							&& !__clashes_ligand(*pneighbor, curr_conformation, curr_state)) {
							// IMPROVEMENT: clashes_receptor & clashes_ligand might
							// be replaced by energy test ? (maybe not)
							const double nb_ene = __score.non_bonded_energy(pneighbor->get_atoms());
							const double torsion_ene = torsion_energy(curr_state, *pneighbor);
							const double branch_ene = curr_conformation.second
								+ nb_ene 
								+ torsion_ene;
							LinkEnergy next_conformation = curr_conformation;
							next_conformation.first.push_back(pneighbor);
							next_conformation.second = branch_ene;
							dbgmsg("accepting state " << *pneighbor);
							dbgmsg("before insert openset.size() = " << openset.size());
							openset.insert(next_conformation);
							dbgmsg("after insert openset.size() = " << openset.size());
						}
					}
				}
			}
		}
		dbgmsg("ASTAR FINISHED");
		if (min_energy != MAX_ENERGY) {
			dbgmsg("SUCCESS minimum energy conformation at step = " 
				<< iter << " = " << endl << min_conformation);
		} else {
			dbgmsg("FAILED to connect start conformation : " << start_conformation);
			throw ConnectionError("die : could not connect this conformation",
				failed_state_pairs);
		}
		Conf mc;
		for (auto &pstate : min_conformation.first) mc.first.push_back(*pstate);
		mc.second = min_conformation.second;
		return mc;
	}
	
	//~ Array2d<bool> Linker::__find_compatible_state_pairs(const State::Vec &states) {
		//~ /* Find all pairs of compatible states at correct distances for 
		 //~ * the multi-seed molecules
		 //~ * 
		 //~ */
		//~ int sz = states.size();
		//~ // reserve memory for sz * sz adjacency matrix; sz...sum of all docked states
		//~ Array2d<bool> conn(sz, sz);
//~ 
		//~ for (int i = 0; i < states.size(); ++i) {
			//~ State &state1 = *states[i];
			//~ const Segment &segment1 = state1.get_segment();
			//~ for (int j = i + 1; j < states.size(); ++j) {
				//~ State &state2 = *states[j];
				//~ const Segment &segment2 = state2.get_segment();
				//~ if (&segment1 != &segment2) {
					//~ const double max_linker_length = segment1.get_max_linker_length(segment2);
					//~ // const double max_dist = __tol_max_coeff * max_linker_length;
					//~ // const double min_dist = __tol_min_coeff * max_linker_length;
					//~ const double max_dist = max_linker_length + 2.0;
					//~ const double min_dist = max_linker_length - 2.0;
					//~ const Bond &excluded = (segment1.is_adjacent(segment2) 
						//~ ? segment1.get_bond(segment2) : Bond());
//~ #ifndef NDEBUG
					//~ dbgmsg("segment1 = " << segment1 << endl << "segment2 = " << segment2);
					//~ dbgmsg("state1 = " << state1 << endl << "state2 = " << state2);
					//~ if (segment1.is_adjacent(segment2)) dbgmsg("excluded bond is " << excluded);
//~ #endif
					//~ const double dist = __distance(state1, state2);
					//~ if (dist < max_dist && dist > min_dist
						//~ && !state1.clashes(state2, excluded)) {
						//~ dbgmsg("compatible states pair belongs to segments " 
							//~ << segment1.get_seed_id() << " and " 
							//~ << segment2.get_seed_id() << " mll is "
							//~ << max_linker_length
							//~ << " dist is " << dist);
						//~ conn.data[i][j] = conn.data[j][i] = true;
					//~ }
				//~ }
			//~ }
		//~ }
		//~ return conn;
	//~ }
//~ 
	//~ Array2d<bool> Linker::__find_compatible_state_pairs(const State::Vec &states) {
		//~ /* Find all pairs of compatible states at correct distances for 
		 //~ * the multi-seed molecules
		 //~ */
		//~ Array2d<bool> conn(states.size());
		//~ Linker::Poses poses(states);
		//~ for (auto &state : states) {
			//~ State::Vec join_states = poses.get_join_states(state);
			//~ State::Set clashed_states = poses.get_clashed_states(state);
			//~ State::Vec compatible_states = join_states - clashed_states;
			//~ const int i = state.get_id();
			//~ for (auto &state2 : compatible_states) {
				//~ const int j = state2.get_id();
				//~ conn.data[i][j] = conn.data[j][i] = true;
			//~ }
		//~ }
		//~ return conn;
	//~ }
//~ 
	Array2d<bool> Linker::__find_compatible_state_pairs(const Seed::Graph &seed_graph) {
		/* Find all pairs of compatible states at correct distances for 
		 * the multi-seed molecules
		 */
		Array2d<bool> conn(states.size());
		Linker::Poses poses(states);
		for (int u = 0; u < seed_graph.size(); ++u) {
			Segment &segment1 = seed_graph[u].get_segment();
			for (int v = u + 1; v < seed_graph.size(); ++v) {
				Segment &segment2 = seed_graph[v].get_segment();
				JoinAtoms join_atoms(segment1, segment2);
				double max_linker_length = segment1.get_max_linker_length(segment2);
				for (auto &state1 : segment1.get_states()) {
					State::Vec join_states = poses.get_join_states(state, max_linker_length);
					State::Set clashed_states = poses.get_clashed_states(state);
					State::Vec compatible_states = join_states - clashed_states;
					const int i = state.get_id();
					for (auto &state2 : compatible_states) {
						const int j = state2.get_id();
						conn.data[i][j] = conn.data[j][i] = true;
					}
				}
			}
		}
		return conn;
	}

	vector<Linker::LinkEnergy> Linker::__generate_rigid_conformations(const Seed::Graph &seed_graph) {
		Benchmark::reset();
		cout << "Generating rigid conformations of states..." << endl;

		State::Vec states;
		for (auto &seed : seed_graph) 
			for (auto &pstate : seed.get_segment().get_states())
				states.push_back(&*pstate);

		// init adjacency matrix (1 when states are compatible, 0 when not)
		Array2d<bool> conn = __find_compatible_state_pairs(states);

		cout << "find_compatible_state_pairs took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		
		// find all maximum cliques, of size k; k...number of seed segments
		Maxclique m(conn.data, conn.szi);
		vector<vector<int>> qmaxes = m.mcq(seed_graph.size());

		cout << "max clique search took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;

		if (qmaxes.empty())
			throw Error("die : couldn't find any possible conformations for ligand " 
				+ __ligand.name());

		// score conformations by summing up individual fragment energies
		vector<LinkEnergy> possibles_w_energy;
		for (auto &qmax : qmaxes) {
			double energy = 0;
			State::Vec conf;
			for (auto &i : qmax) {
				conf.push_back(states[i]);
				energy += states[i]->get_energy();
			}
			possibles_w_energy.push_back({conf, energy});
		}
		
		// sort conformations according to their total energy
		sort(possibles_w_energy.begin(), possibles_w_energy.end(), 
			[] (const LinkEnergy &i, const LinkEnergy &j) { 
			return i.second < j.second;	});
			
		if (__max_possible_conf != -1 && possibles_w_energy.size() > __max_possible_conf) {
			possibles_w_energy.resize(__max_possible_conf);
			dbgmsg("number of possible conformations > max_possible_conf, "
				<< "resizing to= " << __max_possible_conf << " conformations");
		}

		//~ cout << possibles_w_energy << endl;

		dbgmsg("RIGID CONFORMATIONS FOR LIGAND " << __ligand.name() 
			<< " : " << endl << possibles_w_energy);

		cout << "Generated " << possibles_w_energy.size() 
			<< " possible conformations for ligand " << __ligand.name()
			<< ", which took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return possibles_w_energy;
	}
	
	Molecules Linker::connect() {
		Benchmark::reset();
		cout << "Starting connection of seeds for ligand " << __ligand.name() << endl;
		
		Segment::Graph segment_graph = Segment::create_graph(__ligand);
		if (!segment_graph.find_cycles_connected_graph().empty()) {
			throw Error("die : cyclic molecules are currently not supported");
		}
		dbgmsg("segment graph for ligand " << __ligand.name() << " = " << segment_graph);
		__create_states(segment_graph, __top_seeds);

		const Segment::Paths paths = __find_paths(segment_graph);
		__init_max_linker_length(paths);
		__set_branching_rules(paths);
		
		const Seed::Graph seed_graph = Seed::create_graph(segment_graph, paths);
		dbgmsg("seed graph for ligand " << __ligand.name() << " = " << seed_graph);
		
		const vector<LinkEnergy> possible_states = __generate_rigid_conformations(seed_graph);
		Conformations good_conformations = __connect(segment_graph.size(), possible_states);
		
		cout << "Connection of seeds for ligand " << __ligand.name() 
			<< " resulted in " << good_conformations.size() 
			<< " conformations and took " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		return __reconstruct(good_conformations);
	}

	bool Linker::__has_blacklisted(const State::Vec &conformation, 
		const set<State::ConstPair> &blacklist) {

		for (int i = 0; i < conformation.size(); ++i) {
			for (int j = i + 1; j < conformation.size(); ++j) {
				if (blacklist.count({conformation[i], conformation[j]})
					|| blacklist.count({conformation[j], conformation[i]}))
					return true;
			}
		}
		return false;
	}

	Linker::Conformations Linker::__connect(const int segment_graph_size, 
		const vector<LinkEnergy> &possibles) {
		
		Conformations good_conformations;
		set<State::ConstPair> blacklist;
		for (auto &le : possibles) {
			auto &conformation = le.first;
			auto &energy = le.second;
			if (!__has_blacklisted(conformation, blacklist)) {
				LinkEnergy grow_link_energy(conformation, energy);
				dbgmsg("connecting ligand " << __ligand.name() 
					<< endl << "possible starting conformation is " << endl
					<< grow_link_energy);
				try {
					good_conformations.push_back(__a_star(segment_graph_size, 
						grow_link_energy, __link_iter));
					dbgmsg("complete linked molecule" << endl 
						<< good_conformations.back());
				} catch (ConnectionError& e) {
					dbgmsg("exception : " << e.what());
					for (auto &failed_pair : e.get_failed_state_pairs())
						blacklist.insert(failed_pair);
				}
			}
		}
		return good_conformations;
	}
	
	Segment::Paths Linker::__find_paths(const Segment::Graph &segment_graph) {
		/*
		 * Find ALL paths between ALL seed segments (even non-adjacent) 
		 * and seeds and leafs
		 */
		Segment::Paths paths;
		dbgmsg("find all paths in a graph");
		Segment::Set seeds, seeds_and_leafs;
		for (auto &seg : segment_graph) {
			if (seg.is_seed())
				seeds.insert(&seg);
			if (seg.is_seed() || seg.is_leaf())
				seeds_and_leafs.insert(&seg);
		}
		set<Segment::ConstPair> visited;
		for (auto &pseg1 : seeds) {
			for (auto &pseg2 : seeds_and_leafs) {
				if (pseg1 != pseg2 && !visited.count({pseg1, pseg2})) {
					dbgmsg("finding path between " << *pseg1 << " and "
						<< *pseg2);
					visited.insert({pseg1, pseg2});
					visited.insert({pseg2, pseg1});
					Segment::Graph::Path path = Glib::find_path(*pseg1, *pseg2);
					paths.insert({{pseg1, pseg2}, path });
					dbgmsg("path between first segment " << *pseg1
						<< " and second segment " << *pseg2);
				}
			}
		}
		return paths;
	}
	
	void Linker::__set_branching_rules(const Segment::Paths &paths) {
		for (auto &kv : paths) {
			const Segment::Graph::Path &path = kv.second;
#ifndef NDEBUG
			dbgmsg("valid_path = ");
			for (auto &seg : path) dbgmsg(*seg);
#endif
			Segment &start = *path.back();
			Segment &goal = *path.front();
			Segment &start_next = **(path.end() - 2);
			Segment &goal_next = **(path.begin() + 1);
			for (auto it = path.begin(); it != path.end(); ++it) {
				Segment &current = **it;
				if (&current != &goal && goal.is_seed()) {
					Segment &next = **(it - 1);
					current.set_next(goal, next);
					goal.set_next(current, goal_next);
					current.set_adjacent_seed_segments(goal);
					dbgmsg("current = " << current << " goal.is_seed() = " 
						<< boolalpha << goal.is_seed() << " next = " << next);
				}
				if (&current != &start && start.is_seed()) {
					Segment &prev = **(it + 1);
					current.set_next(start, prev);
					start.set_next(current, start_next);
					current.set_adjacent_seed_segments(start);
					dbgmsg("current = " << current << " start.is_seed() = " 
						<< boolalpha << start.is_seed() << " prev = " << prev);
				}
			}
		}
	}
	
	Molecules Linker::__reconstruct(const Conformations &conformations) {
		Benchmark::reset();
		cout << "Reconstructing docked ligands..." << endl;
		Molecules mols;
		int conf_number = 0;
		for (auto &conformation : conformations) { // write out reconstructed molecules
			for (auto &state : conformation.first) { // convert ligand to new coordinates
				const Segment &segment = state.get_segment();
				AtomSet overlap;
				// deal with atoms that overlap between states
				for (auto &adjacent : segment) {
					const Bond &b = segment.get_bond(adjacent);
					overlap.insert(&b.atom2());
				}
				for (auto &kv : state.get_atoms()) {
					Atom &atom = const_cast<Atom&>(*kv.first); // ugly, correct this
					if (!overlap.count(&atom)) {
						const Geom3D::Coordinate &crd = kv.second;
						atom.set_crd(crd);
					}
				}
			}
			Molecule &added = mols.add(new Molecule(__ligand));
			added.set_name(added.name() + "_" + help::to_string(++conf_number));
		}
		cout << "Reconstruction of molecules took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return mols;
	}
};
