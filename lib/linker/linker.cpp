#include "linker.hpp"
#include "geom3d/quaternion.hpp"
#include "score/score.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
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
	Linker::SegGraph Linker::__create_segment_graph(const Molecule &molecule) {
		dbgmsg("Create segment graph ...");
		const Model &model = molecule.first().first();
		vector<unique_ptr<Segment>> vertices;
		for (auto &kv : model.get_rigid()) { // make vertices (segments) of a graph
			const string &nm = kv.first;
			const set<AtomSet> &atom_sets = kv.second;
			dbgmsg(nm);
			for (auto &atoms : atom_sets) {
				vertices.push_back(unique_ptr<Segment>(new Segment(atoms)));
			}
		}
		// connect segments
		for (int i = 0; i < vertices.size(); ++i) {
			Segment &s1 = *vertices[i];
			for (int j = i + 1; j < vertices.size(); ++j) {
				Segment &s2 = *vertices[j];
				auto inter = Glib::intersection(s1.get_atoms(), s2.get_atoms());
				if (inter.size() == 2) {
					s1.add(&s2);
					s2.add(&s1);
					auto &atom1 = **inter.begin();
					auto &atom2 = **inter.rbegin();
					int num_bonds = 0; 
					//~ for (auto &adj : atom1) if (s1.has_atom(adj)) { num_bonds++; }
					//~ for (auto &bond : atom1) {
						//~ auto &adj = bond.second_atom();
					for (auto &adj : atom1) {
						if (s1.has_atom(adj))
							num_bonds++;
					}
					if (num_bonds == 1) {
						dbgmsg("atom " << atom1 << " belongs to segment " << s2);
						dbgmsg("atom " << atom2 << " belongs to segment " << s1);
						s1.set_bond(s2, Bond(&atom2, &atom1));
						s2.set_bond(s1, Bond(&atom1, &atom2));
					} else {
						dbgmsg("atom " << atom1 << " belongs to segment " << s1);
						dbgmsg("atom " << atom2 << " belongs to segment " << s2);
						s1.set_bond(s2, Bond(&atom1, &atom2));
						s2.set_bond(s1, Bond(&atom2, &atom1));
					}
				}
			}
		}
		// find out which rigid segments are seeds
		for (auto &kv : model.get_seeds()) { 
			const string &nm = kv.first;
			//~ dbgmsg("nm = " << nm << " seed_atoms = " << endl << kv.second
			dbgmsg("nm = " << nm);
			for (auto &seed_atoms : kv.second) {
				dbgmsg("seed atoms = " << seed_atoms);
				// loop over segments and find the largest one
				int max_sz = 0;
				Segment *largest;
				for (auto &psegment : vertices) {
					dbgmsg("segment atoms = " << psegment->get_atoms());
					auto inter = Glib::intersection(seed_atoms, psegment->get_atoms());
					if (inter.size() > max_sz) {
						max_sz = inter.size();
						largest = &*psegment;
					}
				}
				largest->set_name(nm);
				largest->set_seed_atoms(seed_atoms);
				dbgmsg("largest_segment_in_seed is " << *largest);
			}
		}
		return SegGraph(std::move(vertices), true, false);
	}
	bool Linker::__link_adjacent(const Linker::SegGraph::Path &path) {
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
	Linker::SeedGraph Linker::__create_seed_graph(const SegGraph &segment_graph,
		const Paths &paths) {
		dbgmsg("Create seed graph ...");
		vector<unique_ptr<Seed>> vertices;
		for (auto &s : segment_graph) // add vertices to seed graph
			if (s.is_seed()) 
				vertices.push_back(unique_ptr<Seed>(new Seed(s)));
		for (int i = 0; i < vertices.size(); ++i) { // ... edges ...
			for (int j = i + 1; j < vertices.size(); ++j) {
				const Segment &seg_i = vertices[i]->get_segment();
				const Segment &seg_j = vertices[j]->get_segment();
				dbgmsg("i = " << i << " j = " << j << " "
					<< seg_i << " " << seg_j << " path exists = "
					<< boolalpha << (paths.count({&seg_i, &seg_j})
					|| paths.count({&seg_j, &seg_i})));
				if (paths.count({&seg_i, &seg_j}) || paths.count({&seg_j, &seg_i})) {
#ifndef NDEBUG
					SegGraph::Path path = paths.count({&seg_i, &seg_j}) ? 
						paths.at({&seg_i, &seg_j}) : paths.at({&seg_j, &seg_i});
					dbgmsg("path = ");
					for (auto it = path.begin(); it != path.end(); ++it)
						dbgmsg(**it);
#endif
					Seed &seed1 = *vertices[i];
					Seed &seed2 = *vertices[j];
					seed1.add(&seed2);
					seed2.add(&seed1);
				}
			}
		}
		return SeedGraph(std::move(vertices), true, false);
	}
	void Linker::__init_max_linker_length(const Paths &paths) {
		for (auto &kv : paths) {
			auto &seg_pair = kv.first;
			SegGraph::Path path(kv.second.begin(), kv.second.end());
			__compute_max_linker_length(path);
			reverse(path.begin(), path.end());
			__compute_max_linker_length(path);
		}
	}
	void Linker::__compute_max_linker_length(SegGraph::Path &path) {
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
	void Linker::__create_states(const SegGraph &segment_graph, const NRset &top_seeds) {
		/* each state is a mapping of docked seed atoms to ligands's 
		 * segment atoms
		 */
		for (auto &seed_mols : top_seeds)
		for (auto &seed_molecule : seed_mols) {
			const string nm = seed_molecule.name();
			dbgmsg("seed name = " << nm);
			/* loop over segments with this name (only those that have 
			 * links to outside), seeds coordinates will not change
			 */
			for (auto &segment : segment_graph) {
				dbgmsg("segment name = " << segment.get_name());
				if (segment.get_name() == nm) {
					// create a graph out of the seed "nm" of "to be" ligand molecule
					MolGraph gs = create_graph(seed_molecule.first().first().get_atoms());
					dbgmsg("gs = " << endl << gs);
					MolGraph gd = create_graph(segment.get_seed_atoms());
					dbgmsg("gd = " << endl << gd);
					// map seed molecules atom numbers to ligand molecule
					MolGraph::Matches m = gd.match(gs);
					dbgmsg("__create_states : m.size() = " << m.size());
#ifndef NDEBUG
					MolGraph::Matches m2 = gs.match(gd);
					dbgmsg("__create_states : m2.size() = " << m2.size());
#endif
					// for each mapping of seed to seed ...
					for (auto &mv : m) {
						// create a state
						auto &vertices1 = mv.first;
						auto &vertices2 = mv.second;
						dbgmsg("new match");
						AtomToCrd atom_crd;
						for (int i = 0; i < vertices2.size(); ++i) {
							Atom &v1 = gd[vertices1[i]];
							Atom &v2 = gs[vertices2[i]];
							//~ Bond &v1 = gd[vertices1[i]];
							//~ Bond &v2 = gs[vertices2[i]];
							if (segment.has_atom(v1)) { // keep only segment atoms
								atom_crd[&v1] = Geom3D::Coordinate(v2.crd());
								dbgmsg("adding matched vertex pair " << vertices1[i] 
									<< "," << vertices2[i] << " new coordinates of atom " 
									<< v1.atom_number() << " are " 
									<< Geom3D::Coordinate(v2.crd()));
							}
#ifndef NDEBUG							
							else {
								dbgmsg("not adding matched vertex pair " << vertices1[i] 
									<< "," << vertices2[i] << " new coordinates of atom " 
									<< v1.atom_number() << " are " 
									<< Geom3D::Coordinate(v2.crd()));
							}
#endif
						}
							//~ if (segment.has_atom(v1.first_atom()) 
								//~ && segment.has_atom(v1.second_atom())) { // keep only segment atoms
								//~ atom_crd[&v1.first_atom()] = Geom3D::Coordinate(v2.first_atom().crd());
								//~ atom_crd[&v1.second_atom()] = Geom3D::Coordinate(v2.second_atom().crd());
								//~ dbgmsg("adding matched vertex pair " << vertices1[i] 
									//~ << "," << vertices2[i] << " new coordinates of atom " 
									//~ << v1.first_atom().atom_number() << " are " 
									//~ << Geom3D::Coordinate(v2.first_atom().crd()) 
									//~ << " and new coordinates of atom "
									//~ << v1.second_atom().atom_number() << " are " 
									//~ << Geom3D::Coordinate(v2.second_atom().crd()));
							//~ }
//~ #ifndef NDEBUG							
							//~ else {
								//~ dbgmsg("not adding matched vertex pair " << vertices1[i] 
									//~ << "," << vertices2[i] << " new coordinates of atom " 
									//~ << v1.first_atom().atom_number() << " are " 
									//~ << Geom3D::Coordinate(v2.first_atom().crd()) 
									//~ << " and new coordinates of atom "
									//~ << v1.second_atom().atom_number() << " are " 
									//~ << Geom3D::Coordinate(v2.second_atom().crd()));
							//~ }
//~ #endif
						//~ }
						// ONLY COPY SEGMENT COORDS NOT SEED
						segment.add_state(unique_ptr<State>(new State(segment, atom_crd, 
							__score.non_bonded_energy(atom_crd))));
					}
				}
			}
		}
		/* unseed those seed segments that don't have states - small 
		 * segments within seed segments - these have no docked coordinates 
		 * (no states)
		 */
		for (auto &segment : segment_graph) {
			if (segment.is_seed() && segment.get_states().empty()) {
				dbgmsg("no representative docked seeds were found for segment name = " 
					<< segment << " therefore this segment is unseeded...");
				segment.set_name("");
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
	
	StateVec Linker::__compute_neighbors(const State &curr_state, Segment &next,
		vector<unique_ptr<State>> &states) {
		dbgmsg("in compute_neighbors for state = " << curr_state);
		StateVec l;
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

		SegmentSet done_segments, free_seeds;
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
		set<ConstStatePair>	failed_state_pairs;
		LinkEnergy min_conformation;
		double min_energy = MAX_ENERGY;
		PriorityQueue openset; // openset has conformations sorted from lowest to highest energy
		openset.insert(LinkEnergy{StateVec{&*states[0]}, states[0]->get_energy()});
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
					auto ret2 = (adj_is_seed ? StateVec{adj_is_seed} : 
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
	
	map<State*, StateSet> Linker::__find_compatible_state_pairs(const SeedGraph &seed_graph) {
		/* Find all pairs of compatible states at correct distances for 
		 * the multi-seed molecules
		 * 
		 */
		set<const Seed*> visited;
		map<State*, StateSet> pos;
		for (auto &seed1 : seed_graph) {
			const Segment &segment1 = seed1.get_segment();
			visited.insert(&seed1);
			for (auto &pstate1 : segment1.get_states()) {
				State &state1 = *pstate1;
				for (auto &seed2 : seed1) {
					if (!visited.count(&seed2)) {
						const Segment &segment2 = seed2.get_segment();
						const double max_linker_length = segment1.get_max_linker_length(segment2);
						const double max_dist = __tol_max_coeff * max_linker_length;
						const double min_dist = __tol_min_coeff * max_linker_length;
						const Bond &excluded = (segment1.is_adjacent(segment2) 
							? segment1.get_bond(segment2) : Bond());
						dbgmsg("segment1 = " << segment1 << endl
							<< "segment2 = " << segment2);
#ifndef NDEBUG
						if (segment1.is_adjacent(segment2))
							dbgmsg("excluded bond is " << excluded);
#endif
						for (auto &pstate2 : segment2.get_states()) {
							State &state2 = *pstate2;
							const double dist = __distance(state1, state2);
							if (dist < max_dist && dist > min_dist
								&& !state1.clashes(state2, excluded)) {
								dbgmsg("compatible states pair belongs to segments " 
									<< segment1.get_name() << " and " 
									<< segment2.get_name() << " mll is "
									<< max_linker_length
									<< " dist is " << dist);
								pos[&state1].insert(&state2);
							}
						}
					}
				}
			}
		}
		return pos;
	}

	vector<vector<StateVec>> Linker::__grow_possibles(const map<State*, StateSet> &pos) {
		/* Create all possible combinations of states of size 2, 3 and 4
		 * and in the case of 3 and 4 seed conformations check for cliques
		 * and that all combinations of states are possible; check also
		 * for atom-clashes between states that aren't seed-adjacent
		 */
		vector<vector<StateVec>> possibles(3);
		auto &possibles2 = possibles[0];
		auto &possibles3 = possibles[1];
		auto &possibles4 = possibles[2];
		dbgmsg("in grow_possibles");
		for (auto &kv : pos) {
			State *pstate1 = kv.first;
			if (!pos.count(pstate1)) continue;
			const Segment &segment1 = pstate1->get_segment();
			for (auto &pstate2 : pos.at(pstate1)) {
				if (!pos.count(pstate2)) continue;
				const Segment &segment2 = pstate2->get_segment();
				possibles2.push_back({pstate1, pstate2});
				for (auto &pstate3 : pos.at(pstate2)) {
					if (!pos.count(pstate3)) continue;
					const Segment &segment3 = pstate3->get_segment();
					 // check : cliques of 3 need to be closed (all pairs
					 // of states must be possible)
					if (segment1.is_seed_adjacent(segment3)
						&& !(pos.at(pstate1).count(pstate3) || pos.at(pstate3).count(pstate1)))
						continue; // fail : state1 and state3 are not possible
					const Bond &excluded13 = (segment1.is_adjacent(segment3) 
						? segment1.get_bond(segment3) : Bond());
					if (pstate1->clashes(*pstate3, excluded13)) continue;
					possibles3.push_back({pstate1, pstate2, pstate3});
					for (auto &pstate4 : pos.at(pstate3)) {
						if (!pos.count(pstate4)) continue;
						const Segment &segment4 = pstate4->get_segment();
						 // check : cliques of 4 need to be closed (again all
						 // pairs of states must be possible
						if (segment1.is_seed_adjacent(segment4)
							&& !(pos.at(pstate1).count(pstate4) || pos.at(pstate4).count(pstate1)))
							continue; // fail : state1 and state4 are not possible
						if (segment2.is_seed_adjacent(segment4)
							&& !(pos.at(pstate2).count(pstate4) || pos.at(pstate4).count(pstate2)))
							continue; // fail : state2 and state4 are not possible
						const Bond &excluded14 = (segment1.is_adjacent(segment4) 
							? segment1.get_bond(segment4) : Bond());
						const Bond &excluded24 = (segment2.is_adjacent(segment4) 
							? segment2.get_bond(segment4) : Bond());
						if (pstate1->clashes(*pstate4, excluded14) 
							|| pstate2->clashes(*pstate4, excluded24)) continue;
						possibles4.push_back({pstate1, pstate2, pstate3, pstate4});
						dbgmsg("just added four-state : " << endl << possibles4.back());
					}
				}
			}
		}
		dbgmsg("out of grow_possibles");
		return possibles;
	}

	vector<Linker::LinkEnergy> Linker::__find_possible_states(const SeedGraph &seed_graph) {
		Benchmark::reset();
		cout << "Finding possible conformations of states..." << endl;
		vector<StateVec> possibles;
		// for one-seed molecule all states are possible
		if (seed_graph.size() == 1) {
			for (auto &pstate : seed_graph[0].get_segment().get_states()) {
				possibles.push_back({&*pstate});
			}
		// for multi-seed molecules
		} else {
			auto pos = __find_compatible_state_pairs(seed_graph);
			auto all_possibles = __grow_possibles(pos);
			
			cout << "Report for possible conformations of ligand " 
				<< __ligand.name() << " : " << endl 
				<< "# of 2-state conformations : " << all_possibles[0].size() << endl 
				<< "# of 3-state conformations : " << all_possibles[1].size() << endl 
				<< "# of 4-state conformations : " << all_possibles[2].size() << endl;
#ifndef NDEBUG
			for (auto &conf : all_possibles[2]) {
				dbgmsg("conformation : ");
				for (auto &pstate : conf) {
					dbgmsg(pstate->pdb());
				}
			}
#endif
			possibles.insert(possibles.end(), all_possibles[2].begin(), all_possibles[2].end());
			possibles.insert(possibles.end(), all_possibles[1].begin(), all_possibles[1].end());
			possibles.insert(possibles.end(), all_possibles[0].begin(), all_possibles[0].end());
		}
		// sort possible conformations according to their energies 
		// (lowest energies are first)
		vector<LinkEnergy> possibles_w_energy;
		for (auto &conformation : possibles) {
			double energy = 0;
			for (auto &pstate : conformation) energy += pstate->get_energy();
			possibles_w_energy.push_back({conformation, energy});
		}
		
		sort(possibles_w_energy.begin(), possibles_w_energy.end(), 
			[] (const LinkEnergy &i, const LinkEnergy &j) { 
			return i.second < j.second;	});
#ifndef NDEBUG
		for (auto &le : possibles_w_energy) {
			auto &conformation = le.first;
			auto &energy = le.second;
			dbgmsg("++++++++ this is one possible conformation for ligand "
				<< __ligand.name() << " of size " << conformation.size() 
				<< " and energy of " << setprecision(4) << fixed << energy);
			for (auto &s : conformation) dbgmsg("    " << *s);
		}			
#endif
		dbgmsg("possibles size = " << possibles_w_energy.size());
		if (__max_possible_conf != -1 && possibles_w_energy.size() > __max_possible_conf) {
			possibles_w_energy.resize(__max_possible_conf);
			dbgmsg("possible conformations is too large, resizing to= " 
				<< __max_possible_conf << " conformations");
		}
		if (possibles_w_energy.empty())
			throw Error("die : couldn't find any possible conformations for ligand " 
				+ __ligand.name());
		cout << "found " << possibles_w_energy.size() 
			<< " possible conformations for ligand " << __ligand.name()
			<< ", which took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return possibles_w_energy;
	}
	
	Molecules Linker::connect() {
		Benchmark::reset();
		cout << "Starting connection of seeds for ligand " << __ligand.name() << endl;
		
		SegGraph segment_graph = __create_segment_graph(__ligand);
		if (!segment_graph.find_cycles_connected_graph().empty()) {
			throw Error("die : cyclic molecules are currently not supported");
		}
		dbgmsg("segment graph for ligand " << __ligand.name() << " = " << segment_graph);
		__create_states(segment_graph, __top_seeds);
		const Paths paths = __find_paths(segment_graph);
		__init_max_linker_length(paths);
		__set_branching_rules(paths);
		const SeedGraph seed_graph = __create_seed_graph(segment_graph, paths);
		dbgmsg("seed graph for ligand " << __ligand.name() << " = " << seed_graph);
		
		Conformations good_conformations = __connect(segment_graph.size(), 
			__find_possible_states(seed_graph));
		
		cout << "Connection of seeds for ligand " << __ligand.name() 
			<< " resulted in " << good_conformations.size() 
			<< " conformations and took " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		return __reconstruct(good_conformations);
	}

	bool Linker::__has_blacklisted(const StateVec &conformation, 
		const set<ConstStatePair> &blacklist) {

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
		set<ConstStatePair> blacklist;
		for (auto &le : possibles) {
			auto &conformation = le.first;
			auto &energy = le.second;
			if (!__has_blacklisted(conformation, blacklist)) {
				//~ LinkEnergy test_conf;
				//~ for (auto &pstate : conformation) {
					//~ if (pstate->get_segment().get_name() == "2")
							//~ test_conf.first.push_back(pstate);
							//~ test_conf.second = pstate->get_energy();
				//~ }
				//~ if (!test_conf.first.empty()) {
					//~ LinkEnergy grow_link_energy(test_conf.first, test_conf.second);
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
					//~ break;
				//~ }
			}
		}
		return good_conformations;
	}
	Linker::Paths Linker::__find_paths(const SegGraph &segment_graph) {
		Paths paths;
		dbgmsg("find all paths in a graph");
		set<Segment*> seeds, seeds_and_leafs;
		for (auto &seg : segment_graph) {
			if (seg.is_seed())
				seeds.insert(&seg);
			if (seg.is_seed() || seg.is_leaf())
				seeds_and_leafs.insert(&seg);
		}
		set<ConstSegPair> visited;
		for (auto &pseg1 : seeds) {
			for (auto &pseg2 : seeds_and_leafs) {
				if (pseg1 != pseg2 && !visited.count({pseg1, pseg2})) {
					dbgmsg("finding path between " << *pseg1 << " and "
						<< *pseg2);
					visited.insert({pseg1, pseg2});
					visited.insert({pseg2, pseg1});
					SegGraph::Path path = Glib::find_path(*pseg1, *pseg2);
					if (__link_adjacent(path)) {
						paths.insert({{pseg1, pseg2}, path });
						dbgmsg("path between first segment " << *pseg1
							<< " and second segment " << *pseg2);
					}
				}
			}
		}
		return paths;
	}
	void Linker::__set_branching_rules(const Paths &paths) {
		for (auto &kv : paths) {
			const SegGraph::Path &path = kv.second;
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
					//~ overlap.insert(&b.second_atom());
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
