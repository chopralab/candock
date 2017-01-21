#include "linker.hpp"
#include "poses.hpp"
#include "geom3d/quaternion.hpp"
#include "score/score.hpp"
#include "pdbreader/nrset.hpp"
#include "pdbreader/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "graph/mcqd.hpp"
#include "modeler/modeler.hpp"
#include "geom3d/geom3d.hpp"
#include "cluster/greedy.hpp"
#include <queue>
#include <iostream>

#include "cuda_linker.h"

using namespace std;

namespace Linker {

	/**
	 * Iterative minimization starts from each seed in the docked molecule. 
	 * It takes top N percent of each seed's docked conformations, and then
	 * reconstructs the whole molecule from that seed.
	 * 
	 */

	DockedConformation Linker::CUDA_IterativeLinker::__a_star(const int segment_graph_size, 
		const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter) {
		
        cout << "Segment Graph Size " << segment_graph_size << endl;
        cout << "Iter " << iter << endl;
        


		for (auto &pstate : start_conformation.get_states())
			states.push_back(unique_ptr<State>(new State(*pstate)));
		if (start_conformation.empty())
			throw Error ("die : at least one docked anchor state is required for linking");
		SegStateMap docked_seeds;
		for (auto &pstate : states) 
			docked_seeds.insert({&pstate->get_segment(), &*pstate});

        //Move data to GPU
        cuda_linker cuda;
        cuda.setup();

		set<State::ConstPair> failed_state_pairs;
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
	
	/**
	 * For iterative create top percent states of each seed. Here, partial
	 * contains only one state (or any desired number up to seed_graph.size())
	 * 
	 */

	Partial::Vec Linker::CUDA_IterativeLinker::__generate_rigid_conformations(const Seed::Graph &seed_graph) {
		Benchmark bench;
		cout << "Generating rigid conformations of states for " << __ligand.name() << "..." << endl;

		State::Vec states;
		State::Id id = 0;
		for (auto &seed : seed_graph)
			for (auto &pstate : seed.get_segment().get_states()) {
				states.push_back(&*pstate);
				pstate->set_id(id++);
				dbgmsg(*pstate);
			}
				
		//help::memusage("before find compatible state pairs");
	
		// init adjacency matrix (1 when states are compatible, 0 when not)
		Array2d<bool> conn = __find_compatible_state_pairs(seed_graph, states.size());

		dbgmsg("dimensions of conn = " << conn.get_szi() << " " << conn.get_szj());
		dbgmsg("conn = " << conn);

		//help::memusage("before max.clique.search");
		
		cout << "find_compatible_state_pairs took " << bench.seconds_from_start() 
			<< " wallclock seconds for " << __ligand.name() << endl;
		
		Partial::Vec possibles_w_energy;
		{
			// find all maximum cliques with the number of seed segments of __max_clique_size
			Maxclique m(conn);
			
			const int mcq_size = __max_clique_size > seed_graph.size() ? seed_graph.size() : __max_clique_size;
			const vector<vector<unsigned short int>> &qmaxes = m.mcq(mcq_size);
	
			//help::memusage("after max.clique.search");
	
			cout << "found " << qmaxes.size() << " maximum cliques, which took " 
				<< bench.seconds_from_start() << " wallclock seconds for " << __ligand.name() << endl;
	
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
				possibles_w_energy.push_back(Partial(conf, energy));
				dbgmsg("partial = " << possibles_w_energy.back());
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
			<< " possible top percent docked seeds that will serve as starting points for reconstruction of ligand " << __ligand.name()
			<< ", which took " << bench.seconds_from_start() 
			<< " wallclock seconds" << endl;
		return clustered_possibles_w_energy;
	}
	

	DockedConformation Linker::CUDA_IterativeLinker::__reconstruct(const Partial &conformation) {
		Benchmark bench;
		cout << "Reconstructing docked ligands for " << __ligand.name() << "..." << endl;
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
		ligand.set_name(__ligand.name() + "_" + std::to_string(++conf_number));
		
		// receptor
		Molib::Molecule receptor(__receptor);
		const Geom3D::Point::Vec &crds = conformation.get_receptor_crds();
		int i = 0;
		for (auto &patom : receptor.get_atoms()) {
			patom->set_crd(crds[i++]);
		}
		cout << "Reconstruction of molecules took " << bench.seconds_from_start() 
			<< " wallclock seconds for " << __ligand.name() << endl;
		return DockedConformation(ligand, receptor, conformation.get_energy());
	}
};
