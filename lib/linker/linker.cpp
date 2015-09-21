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
#include "modeler/modeler.hpp"
#include "geom3d/geom3d.hpp"
#include "cluster/greedy.hpp"
#include <queue>

using namespace std;

namespace Linker {

	ostream& operator<<(ostream& os, const DockedConformation &conf)	{
		os << "start link ++++++++++++++++++++++++++++++" << endl;
		os << "ENERGY = " << conf.get_energy() << endl;
		os << "LIGAND = " << conf.get_ligand() << endl;
		os << "RECEPTOR = " << conf.get_receptor() << endl;
		os << "end link --------------------------------" << endl;
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
		const Partial &conformation, const State &prev) const {

		const Molib::Bond &excluded = current.get_segment().get_bond(prev.get_segment());
		// clashes between current segment and previous segments
		for (auto &pstate : conformation.get_states()) {
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
		const int start_atom_idx = start_segment
			.get_bond(start_segment.get_next(goal_segment)).idx1();
		const int goal_atom_idx = goal_segment
			.get_bond(goal_segment.get_next(start_segment)).idx1();
		const Geom3D::Coordinate &first_crd = start.get_crd(start_atom_idx);
		const Geom3D::Coordinate &last_crd = goal.get_crd(goal_atom_idx);
		return first_crd.distance(last_crd);
	}
	
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

	string Linker::to_pdb(const Partial &conformation) {
		stringstream os;
		os << "MODEL" << endl;
		os << "REMARK   8 ENERGY " << conformation.get_energy() << endl;
		for (auto &pstate : conformation.get_states())
			os << pstate->pdb();
		os << "ENDMDL" << endl;
		return os.str();
	}

	pair<State*, Segment*> Linker::__find_good_neighbor(
		const Partial &curr_conformation, const SegStateMap &docked_seeds) {

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
	}

	DockedConformation Linker::__a_star_static(const int segment_graph_size, 
		const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter) {

		
		for (auto &pstate : start_conformation.get_states())
			states.push_back(unique_ptr<State>(new State(*pstate)));
		if (start_conformation.empty())
			throw Error ("die : at least one docked anchor state is required for linking");
		SegStateMap docked_seeds;
		for (auto &pstate : states) 
			docked_seeds.insert({&pstate->get_segment(), &*pstate});
		set<State::ConstPair> failed_state_pairs;
		Partial min_conformation(MAX_ENERGY);
		PriorityQueue openset; // openset has conformations sorted from lowest to highest energy
		openset.insert(Partial(State::Vec{&*states[0]}, states[0]->get_energy(),
			Geom3D::Point::Vec()));
			
		while(!openset.empty()) {
			if (--iter < 0) break;
			Partial curr_conformation = *openset.begin();
			openset.erase(openset.begin());
#define MOVIE
#ifdef MOVIE
			DockedConformation docked = __reconstruct(curr_conformation);
			docked.get_receptor().undo_mm_specific();
			inout::output_file(Molib::Molecule::print_complex(docked.get_ligand(), docked.get_receptor(), docked.get_energy()), 
			"movie.pdb", ios_base::app); // output docked molecule conformations
#endif

			dbgmsg("openset.size() = " << openset.size());
			dbgmsg("curr_conformation at step = " 
				<< iter << " = " << endl << curr_conformation
				<< endl << to_pdb(curr_conformation));
			if (curr_conformation.size() == segment_graph_size) {
				dbgmsg("CANDIDATE for minimum energy conformation at step = " 
					<< iter << " = " << endl << curr_conformation);
				if (curr_conformation.get_energy() < min_conformation.get_energy()) {
					min_conformation = curr_conformation;
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
						//~ dbgmsg("clashes_receptor = " << boolalpha 
							//~ << __clashes_receptor(*pneighbor));
						dbgmsg("clashes_ligand = " << boolalpha 
							<< __clashes_ligand(*pneighbor, curr_conformation, curr_state));
						//~ if (!__clashes_receptor(*pneighbor) 
							//~ && !__clashes_ligand(*pneighbor, curr_conformation, curr_state)) {
						if (!__clashes_ligand(*pneighbor, curr_conformation, curr_state)) { // don't check for clashes with receptor
							
							Partial next_conformation(curr_conformation);
							next_conformation.add_state(*pneighbor);

							// compute energy with original (static) receptor structure
							const double energy = __score.non_bonded_energy(__gridrec, next_conformation.get_ligand_atoms(), next_conformation.get_ligand_crds());

							next_conformation.set_energy(energy);

							dbgmsg("accepting state " << *pneighbor);
							openset.insert(next_conformation);
						}
					}
				}
			}
		}
		dbgmsg("ASTAR FINISHED");
		if (min_conformation.get_energy() != MAX_ENERGY) {
			dbgmsg("SUCCESS minimum energy conformation at step = " 
				<< iter << " = " << endl << min_conformation);
		} else {
			dbgmsg("FAILED to connect start conformation : " << start_conformation);
			throw ConnectionError("die : could not connect this conformation",
				failed_state_pairs);
		}
		return __reconstruct(min_conformation);
	}
	
	DockedConformation Linker::__a_star_iterative(const int segment_graph_size, 
		const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter) {
		
		for (auto &pstate : start_conformation.get_states())
			states.push_back(unique_ptr<State>(new State(*pstate)));
		if (start_conformation.empty())
			throw Error ("die : at least one docked anchor state is required for linking");
		SegStateMap docked_seeds;
		for (auto &pstate : states) 
			docked_seeds.insert({&pstate->get_segment(), &*pstate});
		set<State::ConstPair>	failed_state_pairs;
		Partial min_conformation(MAX_ENERGY);
		PriorityQueue openset; // openset has conformations sorted from lowest to highest energy
		openset.insert(Partial(State::Vec{&*states[0]}, states[0]->get_energy(),
			__receptor.get_crds()));
			
		while(!openset.empty()) {
			if (--iter < 0) break;
			Partial curr_conformation = *openset.begin();
			openset.erase(openset.begin());
			dbgmsg("openset.size() = " << openset.size());
			dbgmsg("curr_conformation at step = " 
				<< iter << " = " << endl << curr_conformation
				<< endl << to_pdb(curr_conformation));
			if (curr_conformation.size() == segment_graph_size) {
				dbgmsg("CANDIDATE for minimum energy conformation at step = " 
					<< iter << " = " << endl << curr_conformation);
				if (curr_conformation.get_energy() < min_conformation.get_energy()) {
					min_conformation = curr_conformation;
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
						dbgmsg("clashes_ligand = " << boolalpha 
							<< __clashes_ligand(*pneighbor, curr_conformation, curr_state));
						if (!__clashes_ligand(*pneighbor, curr_conformation, curr_state)) { // don't check for clashes with receptor

							Partial next_conformation(curr_conformation);
							next_conformation.add_state(*pneighbor);

							// prepare for minimization receptor and ligand coordinates
							__modeler.add_crds(__receptor.get_atoms(), next_conformation.get_receptor_crds());
							__modeler.add_crds(next_conformation.get_ligand_atoms(), next_conformation.get_ligand_crds());
							__modeler.init_openmm_positions();
							__modeler.mask(__ligand.get_atoms());
							__modeler.unmask(next_conformation.get_ligand_atoms());
							
#ifndef NDEBUG
							dbgmsg("initial coordinates:");
							for (auto &point : next_conformation.get_ligand_crds()) dbgmsg(point);
#endif

							// minimize ...
							__modeler.minimize_state(__ligand, __receptor, __score);

#ifndef NDEBUG
							dbgmsg("minimized coordinates:");
							for (auto &point : __modeler.get_state(next_conformation.get_ligand_atoms())) dbgmsg(point);
#endif

							// init with minimized coordinates
							Molib::Molecule minimized_receptor(__receptor, __modeler.get_state(__receptor.get_atoms()));
							next_conformation.set_receptor_crds(minimized_receptor.get_crds());
							next_conformation.set_ligand_crds(__modeler.get_state(next_conformation.get_ligand_atoms()));

							// compute energy after minimization
							Molib::Atom::Grid gridrec(minimized_receptor.get_atoms());
							const double energy = __score.non_bonded_energy(gridrec, next_conformation.get_ligand_atoms(), next_conformation.get_ligand_crds());

							next_conformation.set_energy(energy);

							dbgmsg("accepting state " << *pneighbor);
							openset.insert(next_conformation);
						}
					}
				}
			}
		}
		dbgmsg("ASTAR FINISHED");
		if (min_conformation.get_energy() != MAX_ENERGY) {
			dbgmsg("SUCCESS minimum energy conformation at step = " 
				<< iter << " = " << endl << min_conformation);
		} else {
			dbgmsg("FAILED to connect start conformation : " << start_conformation);
			throw ConnectionError("die : could not connect this conformation",
				failed_state_pairs);
		}
		return __reconstruct(min_conformation);
	}
	
	Array2d<bool> Linker::__find_compatible_state_pairs(const Seed::Graph &seed_graph, const int sz) {
		/* Find all pairs of compatible states at correct distances for 
		 * the multi-seed molecules
		 */
		 
		int idx = 0;
		
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
					<< jatoms.second->atom_number() << " excluded = " << excluded);
				for (auto &pstate1 : segment1.get_states()) {
					State &state1 = *pstate1;
					State::Set join_states = poses.get_join_states(state1, segment2, jatoms, 
						max_linker_length, __tol_seed_dist);

					const int i = state1.get_id();
					//~ // dbgmsg("state1(" << i << ") = " << state1);
					for (auto &pstate2 : join_states) {
						//~ // dbgmsg("state2(" << pstate2->get_id() << ") = " << *pstate2);
						if (!state1.clashes(*pstate2, excluded, __clash_coeff)) {
							const int j = pstate2->get_id();
							conn.set(i, j);
							conn.set(j, i);
							//~ dbgmsg("compatible states " << i << " and " << j);
						}
					}
				}
			}
		}

		return conn;
	}

	Partial::Vec Linker::__generate_rigid_conformations(const Seed::Graph &seed_graph) {
		Benchmark::reset();
		cout << "Generating rigid conformations of states..." << endl;

		State::Vec states;
		State::Id id = 0;
		//~ vector<double> scores;
		for (auto &seed : seed_graph)
			for (auto &pstate : seed.get_segment().get_states()) {
				states.push_back(&*pstate);
				//~ scores.push_back(pstate->get_energy());
				pstate->set_id(id++);
				dbgmsg(*pstate);
			}
				
		//help::memusage("before find compatible state pairs");
	
		// init adjacency matrix (1 when states are compatible, 0 when not)
		Array2d<bool> conn = __find_compatible_state_pairs(seed_graph, states.size());

		dbgmsg("dimensions of conn = " << conn.get_szi() << " " << conn.get_szj());

		//help::memusage("before max.clique.search");
		
		cout << "find_compatible_state_pairs took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		
		Partial::Vec possibles_w_energy;
		{
			// find all maximum cliques, of size k; k...number of seed segments
			Maxclique m(conn);
			const vector<vector<unsigned short int>> &qmaxes = m.mcq(seed_graph.size());
	
			//help::memusage("after max.clique.search");
	
			cout << "found " << qmaxes.size() << " maximum cliques, which took " 
				<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
	
			if (qmaxes.empty())
				throw Error("die : couldn't find any possible conformations for ligand " 
					+ __ligand.name());
	
			// score conformations by summing up individual fragment energies
			for (auto &qmax : qmaxes) {
				double energy = 0;
				State::Vec conf;
				for (auto &i : qmax) {
					conf.push_back(states[i]);
					energy += states[i]->get_energy();
				}
				if (energy < __max_allow_energy) {
					possibles_w_energy.push_back(Partial(conf, energy));
				}
			}
	
			//help::memusage("after possibles_w_energy");
		}
		//help::memusage("after forced destructor of qmaxes");
		
		// sort conformations according to their total energy
		Partial::sort(possibles_w_energy);

		if (possibles_w_energy.size() > __max_num_possibles)
			possibles_w_energy.resize(__max_num_possibles);

		//help::memusage("after sort of possibles_w_energy");
		
		dbgmsg("number of possibles_w_energy = " << possibles_w_energy.size());
		// cluster rigid conformations and take only cluster representatives for further linking
		Partial::Vec clustered_possibles_w_energy = Molib::Cluster::greedy(
			possibles_w_energy, __gridrec, __docked_clus_rad);
			
		//help::memusage("after greedy");
			
		dbgmsg("number of clustered_possibles_w_energy = " << clustered_possibles_w_energy.size());
		
		if (__max_possible_conf != -1 && clustered_possibles_w_energy.size() > __max_possible_conf) {
			clustered_possibles_w_energy.resize(__max_possible_conf);
			dbgmsg("number of possible conformations > max_possible_conf, "
				<< "resizing to= " << __max_possible_conf << " conformations");
		}

		dbgmsg("RIGID CONFORMATIONS FOR LIGAND " << __ligand.name() 
			<< " : " << endl << clustered_possibles_w_energy);

		cout << "Generated " << clustered_possibles_w_energy.size() 
			<< " possible conformations for ligand " << __ligand.name()
			<< ", which took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return clustered_possibles_w_energy;
	}
	
	DockedConformation::Vec Linker::connect() {
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
		
		const Partial::Vec possible_states = __generate_rigid_conformations(seed_graph);

		//help::memusage("before connect");

		DockedConformation::Vec docked_confs = __connect(segment_graph.size(), possible_states);

		//help::memusage("after connect");
		
		cout << "Connection of seeds for ligand " << __ligand.name() 
			<< " resulted in " << docked_confs.size() 
			<< " conformations and took " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		return docked_confs;
	}

	Partial::Vec Linker::init_conformations() {
		Benchmark::reset();
		cout << "Starting connection of seeds for ligand " << __ligand.name() << endl;
		
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
	
	DockedConformation Linker::compute_conformation(const Partial &partial) {

		auto &conformation = partial.get_states();
			
		if (!__has_blacklisted(conformation, __blacklist)) {
			try {
				Partial grow_link_energy(conformation, partial.get_energy());
				dbgmsg("connecting ligand " << __ligand.name() 
					<< endl << "possible starting conformation is " << endl
					<< grow_link_energy);
				vector<unique_ptr<State>> states;
				if (__iterative) {
					return __a_star_iterative(__segment_graph.size(), grow_link_energy, states, __link_iter);
				} else {
					return __a_star_static(__segment_graph.size(), grow_link_energy, states, __link_iter);
				}
			} catch (ConnectionError& e) {
				dbgmsg("exception : " << e.what());
				for (auto &failed_pair : e.get_failed_state_pairs())
					__blacklist.insert(failed_pair);
				throw e;
			}
		}
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

	DockedConformation::Vec Linker::__connect(const int segment_graph_size, 
		const Partial::Vec &possibles) {
		
		DockedConformation::Vec good_conformations;
		set<State::ConstPair> blacklist;
		
		for (auto &le : possibles) {
			auto &conformation = le.get_states();
			double energy = le.get_energy();
			
			if (!__has_blacklisted(conformation, blacklist)) {
				try {
					Partial grow_link_energy(conformation, energy);
					dbgmsg("connecting ligand " << __ligand.name() 
						<< endl << "possible starting conformation is " << endl
						<< grow_link_energy);
					vector<unique_ptr<State>> states;
					if (__iterative) {
						good_conformations.push_back(__a_star_iterative(
							segment_graph_size, grow_link_energy, states, __link_iter));
					} else {
						good_conformations.push_back(__a_star_static(
							segment_graph_size, grow_link_energy, states, __link_iter));
					}
					dbgmsg("complete linked molecule" << endl 
						<< good_conformations.back());
					dbgmsg("good_conformations.size() = " << good_conformations.size());
					dbgmsg("good_conformations.capacity() = " << good_conformations.capacity());
					dbgmsg("states.size() = " << states.size());
					dbgmsg("blacklist.size() = " << blacklist.size());
				} catch (ConnectionError& e) {
					dbgmsg("exception : " << e.what());
					for (auto &failed_pair : e.get_failed_state_pairs())
						blacklist.insert(failed_pair);
				}
			}
		}
		return good_conformations;
	}
	
	DockedConformation Linker::__reconstruct(const Partial &conformation) {
		Benchmark::reset();
		cout << "Reconstructing docked ligands..." << endl;
		int conf_number = 0;
		
		// ligand
		for (auto &pstate : conformation.get_states()) { // convert ligand to new coordinates
			State &state = *pstate;
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
		
		Molib::Molecule ligand(__ligand);
		ligand.set_name(__ligand.name() + "_" + help::to_string(++conf_number));
		
		if (__iterative) {
			// receptor
			Molib::Molecule receptor(__receptor);
			const Geom3D::Point::Vec &crds = conformation.get_receptor_crds();
			int i = 0;
			for (auto &patom : receptor.get_atoms()) {
				patom->set_crd(crds[i++]);
			}
			cout << "Reconstruction of molecules took " << Benchmark::seconds_from_start() 
				<< " wallclock seconds" << endl;
			return DockedConformation(ligand, receptor, conformation.get_energy());
		}
		cout << "Reconstruction of molecules took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return DockedConformation(ligand, __receptor, conformation.get_energy());
	}
};
