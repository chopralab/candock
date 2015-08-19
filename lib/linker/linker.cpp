#include "linker.hpp"
#include "poses.hpp"
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

namespace Linker {
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
			//~ os << *pstate << endl;
			os << pstate->pdb() << endl;
		os << "end link --------------------------------" << endl;
		return os;
	}

	ostream& operator<<(ostream& os, const vector<Linker::LinkEnergy> &vec_le)	{
		for (auto &le : vec_le) {
			os << le << endl;
		}			
		return os;
	}

	void Linker::__create_states(const Segment::Graph &segment_graph, const Molib::NRset &top_seeds) {
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
						for (int i = 0; i < vertices2.size(); ++i) {
							crds[vertices1[i]] = seed_atoms.at(vertices2[i])->crd();
							dbgmsg("adding matched vertex pair " << vertices1[i] 
								<< "," << vertices2[i] << " new coordinates of atom " 
								<< segment.get_atom(vertices1[i]).atom_number() << " are " 
								<< crds[vertices1[i]]);
						}
						const double energy = stod(seed_molecule.name());
						dbgmsg("adding docked state " << segment.get_seed_id() << " with energy of " << energy);
						// ONLY COPY SEGMENT COORDS NOT SEED
						segment.add_state(unique_ptr<State>(new State(segment, crds, energy)));
					}
				}
			}
		}
	}
	
	static double torsion_energy(const State &first, const State &second) { return 0.0; }

	bool Linker::__clashes_receptor(const State &current) const {
		for (int i = 0; i < current.get_crds().size(); ++i) {
			const Molib::Atom &a = current.get_segment().get_atom(i); 
			const Geom3D::Coordinate &c = current.get_crd(i); 
			dbgmsg("in clashes_receptor test coordinate = " << c);
			if (__gridrec.clashes(Molib::Atom(c, a.idatm_type()), __clash_coeff)) return true;
		}
		return false;
	}
	
	bool Linker::__clashes_ligand(const State &current, 
		const LinkEnergy &conformation, const State &prev) const {

		const Molib::Bond &excluded = current.get_segment().get_bond(prev.get_segment());
		// clashes between current segment and previous segments
		for (auto &pstate : conformation.first) {
			dbgmsg("clashes_ligand test for state " << current 
				<< " and state " << *pstate);
			if (current.clashes(*pstate, excluded, __clash_coeff)) 
				return true;
		}
		return false;
	}

	double Linker::__distance(const State &start, const State &goal) const {
		const Segment &start_segment = start.get_segment();
		const Segment &goal_segment = goal.get_segment();
		//~ const Molib::Atom &first_atom = start_segment
			//~ .get_bond(start_segment.get_next(goal_segment)).atom1();
		//~ const Molib::Atom &last_atom = goal_segment
			//~ .get_bond(goal_segment.get_next(start_segment)).atom1();
		const int start_atom_idx = start_segment
			.get_bond(start_segment.get_next(goal_segment)).idx1();
		const int goal_atom_idx = goal_segment
			.get_bond(goal_segment.get_next(start_segment)).idx1();
		//~ const Geom3D::Coordinate &first_crd = start.get_atom_crd(first_atom);
		//~ const Geom3D::Coordinate &last_crd = goal.get_atom_crd(last_atom);
		const Geom3D::Coordinate &first_crd = start.get_crd(start_atom_idx);
		const Geom3D::Coordinate &last_crd = goal.get_crd(goal_atom_idx);
		return first_crd.distance(last_crd);
	}
	
	//~ AtomToCrd Linker::__rotate(const Geom3D::Quaternion &q, 
		//~ const Geom3D::Point &p1, const Geom3D::Point &p2, const AtomToCrd &atom_crd) {
		//~ AtomToCrd new_crd;
		//~ for (auto &kv : atom_crd) {	
			//~ new_crd.insert({kv.first, q.rotatedVector(kv.second - p1) + p1}); 
		//~ }
		//~ return new_crd;
	//~ }
	//~ 
	Geom3D::Point::Vec Linker::__rotate(const Geom3D::Quaternion &q, 
		const Geom3D::Point &p1, const Geom3D::Point &p2, const Geom3D::Point::Vec &crds) {

		Geom3D::Point::Vec rotated;
		for (auto &crd : crds) {	
			rotated.push_back(q.rotatedVector(crd - p1) + p1); 
		}
		return rotated;
	}
	
	State::Vec Linker::__compute_neighbors(const State &curr_state, Segment &next,
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
		//~ const Molib::Atom *a1 = &current.adjacent_in_segment(*a2, *a3);	// segments overlap on 1 rotatable bond
		const int idx1 = current.adjacent_in_segment(*a2, *a3);	// segments overlap on 1 rotatable bond
		const Molib::Atom *a1 = &current.get_atom(idx1);
		
		//~ const Molib::Atom *a4 = &next.adjacent_in_segment(*a3, *a2);
		const int idx4 = next.adjacent_in_segment(*a3, *a2);
		const Molib::Atom *a4 = &next.get_atom(idx4);
		
		dbgmsg("a4 = " << *a4 << " coordinates not set yet!");
		Geom3D::Coordinate crd3(curr_state.get_crd(idx3));
		dbgmsg("a3 = " << *a3 << " crd3 = " << curr_state.get_crd(idx3));
		Geom3D::Coordinate crd2(curr_state.get_crd(idx2));
		dbgmsg("a2 = " << *a2 << " crd2 = " << curr_state.get_crd(idx2));
		Geom3D::Coordinate crd1(curr_state.get_crd(idx1));
		dbgmsg("a1 = " << *a1 << " crd1 = " << curr_state.get_crd(idx1));
		states.push_back(unique_ptr<State>(new State(next, 
			__ic.cartesian(*a1, *a2, *a3, crd1, crd2, crd3, next.get_atoms()))));
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
				__rotate(q, crd2, crd3, previous_rotated.get_crds()))));
			dbgmsg("rotated state at angle = " << Geom3D::degrees(angle)
				<< " is " << *states.back());
		}
		State::Vec ret;
		for (int i = ini_sz; i < states.size(); ++i) ret.push_back(&*states[i]);
		return ret;
	}
	
	bool Linker::__check_distances_to_seeds(const State &curr_state, 
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

							/**
							 * FOR MINIMIZATION :
							 * 
							 * openmm.add_state(curr_conformation.get_state());
							 * openmm.add_coords(pneighbor, pneighbor->get_crds());
							 * openmm.minimize_state();
							 * 
							 * double nonbond_energy = score.non_bonded_energy(openmm.get_nonbonded_list(ligand, receptor));
							 * 
							 * 
							 * next_conformation.set_states(curr_conformation.get_states());
							 * next_conformation.add_coords(pneighbor, pneighbor->get_coords());
							 * 
							 * next_conformation.set_state(openmm.get_state());
							 * next_conformation.set_energy(nonbond_energy());
							 * 
							 */
							const double nb_ene = __score.non_bonded_energy(__gridrec, pneighbor->get_segment().get_atoms(), pneighbor->get_crds());
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
	
	Array2d<bool> Linker::__find_compatible_state_pairs(const Seed::Graph &seed_graph, const int sz) {
		/* Find all pairs of compatible states at correct distances for 
		 * the multi-seed molecules
		 */
		Array2d<bool> conn(sz);
		Poses poses(seed_graph);
		for (int u = 0; u < seed_graph.size(); ++u) {
			Segment &segment1 = seed_graph[u].get_segment();
			for (int v = u + 1; v < seed_graph.size(); ++v) {
				Segment &segment2 = seed_graph[v].get_segment();
				const double max_linker_length = segment1.get_max_linker_length(segment2);
				Molib::Atom::Pair jatoms{&segment1.get_bond(segment1.get_next(segment2)).atom1(), 
					&segment2.get_bond(segment2.get_next(segment1)).atom1()};
				const Molib::Bond &excluded = (segment1.is_adjacent(segment2) 
							? segment1.get_bond(segment2) : Molib::Bond());
				dbgmsg("finding compatible state pairs for segments " << segment1 << " and " 
					<< segment2 << " with maximum linker length " << max_linker_length << " and "
					<< " join atom1 = " << jatoms.first->atom_number() << " join atom2 = "
					<< jatoms.second->atom_number());
				for (auto &pstate1 : segment1.get_states()) {
					State &state1 = *pstate1;
					State::Set join_states = poses.get_join_states(state1, segment2, jatoms, 
						max_linker_length, __tol_seed_dist);
					const int i = static_cast<int>(state1.get_id());
					for (auto &pstate2 : join_states) {
						if (!state1.clashes(*pstate2, excluded, __clash_coeff)) {
							const int j = static_cast<int>(pstate2->get_id());
							conn.data[i][j] = conn.data[j][i] = true;
							dbgmsg("compatible states " << i << " and " << j);
						}
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
		Array2d<bool> conn = __find_compatible_state_pairs(seed_graph, states.size());

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
	
	Molib::Molecules Linker::connect() {
		Benchmark::reset();
		cout << "Starting connection of seeds for ligand " << __ligand.name() << endl;
		
		Segment::Graph segment_graph = Segment::create_graph(__ligand);
		if (!segment_graph.find_cycles_connected_graph().empty()) {
			throw Error("die : cyclic molecules are currently not supported");
		}
		dbgmsg("segment graph for ligand " << __ligand.name() << " = " << segment_graph);
		__create_states(segment_graph, __top_seeds);

		const Seed::Graph seed_graph = Seed::create_graph(segment_graph);
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
	
	Molib::Molecules Linker::__reconstruct(const Conformations &conformations) {
		Benchmark::reset();
		cout << "Reconstructing docked ligands..." << endl;
		Molib::Molecules mols;
		int conf_number = 0;
		for (auto &conformation : conformations) { // write out reconstructed molecules
			for (auto &state : conformation.first) { // convert ligand to new coordinates
				const Segment &segment = state.get_segment();
				Molib::Atom::Set overlap;
				// deal with atoms that overlap between states
				for (auto &adjacent : segment) {
					const Molib::Bond &b = segment.get_bond(adjacent);
					overlap.insert(&b.atom2());
				}
				for (int i = 0; i < state.get_segment().get_atoms().size(); ++i) {
					Molib::Atom &atom = const_cast<Molib::Atom&>(state.get_segment().get_atom(i)); // ugly, correct this
					if (!overlap.count(&atom)) {
						const Geom3D::Coordinate &crd = state.get_crd(i);
						atom.set_crd(crd);
					}
				}
			}
			Molib::Molecule &added = mols.add(new Molib::Molecule(__ligand));
			added.set_name(added.name() + "_" + help::to_string(++conf_number));
		}
		cout << "Reconstruction of molecules took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return mols;
	}
};
