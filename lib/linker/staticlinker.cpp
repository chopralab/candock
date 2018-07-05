#include "candock/linker/linker.hpp"
#include "candock/linker/poses.hpp"
#include "candock/geometry/quaternion.hpp"
#include "candock/score/score.hpp"
#include "candock/molib/nrset.hpp"
#include "candock/molib/bond.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/help.hpp"
#include "candock/helper/array2d.hpp"
#include "candock/graph/mcqd.hpp"
#include "candock/modeler/modeler.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/cluster/greedy.hpp"
#include <queue>

using namespace std;

namespace candock {
namespace Linker {

        DockedConformation Linker::StaticLinker::__a_star(const int segment_graph_size, const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter) {

                Benchmark bench;

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
			geometry::Point::Vec()));
			
		while(!openset.empty()) {
			if (--iter < 0) break;
			Partial curr_conformation = *openset.begin();
			openset.erase(openset.begin());
//~ #define MOVIE
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
						dbgmsg("clashes_receptor = " << boolalpha 
							<< __clashes_receptor(*pneighbor));
						dbgmsg("clashes_ligand = " << boolalpha 
							<< __clashes_ligand(*pneighbor, curr_conformation, curr_state));
						//~ if (!__clashes_receptor(*pneighbor) 
							//~ && !__clashes_ligand(*pneighbor, curr_conformation, curr_state)) {
						if (!__clashes_ligand(*pneighbor, curr_conformation, curr_state)) { // don't check for clashes with receptor
							
							Partial next_conformation(curr_conformation);
							next_conformation.add_state(*pneighbor);

							// compute energy with original (static) receptor structure, use expensive non_bonded_energy calculation only for linkers
							const double energy = curr_conformation.get_energy() +
								(pneighbor->get_energy() == 0 ? __score.non_bonded_energy(__gridrec, pneighbor->get_segment().get_atoms(), pneighbor->get_crds()) 
								: pneighbor->get_energy());

							next_conformation.set_energy(energy);
#ifndef NDEBUG
							// TEST if energy determined in a different way is the same
							double ene = 0.0, ene2 = 0.0;
							for (auto &pstate : next_conformation.get_states()) {
								if (pstate->get_energy() == 0) {
									ene += __score.non_bonded_energy(__gridrec, pstate->get_segment().get_atoms(), pstate->get_crds());
								} else {
									ene += pstate->get_energy();
									ene2 += pstate->get_energy();
								}
							}
							dbgmsg("ene = " << ene << " ene2 = " << ene2 << " energy = " << energy);
#endif

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

#ifndef NDEBUG
		DockedConformation mini = __reconstruct(min_conformation);
		dbgmsg("ENERGY OF MINIMIZED LIGAND (reconstructed) = " << __score.non_bonded_energy(__gridrec, mini.get_ligand()));
		dbgmsg("ENERGY OF MINIMIZED LIGAND (calculated) = " << 
			__score.non_bonded_energy(__gridrec, min_conformation.get_ligand_atoms(), min_conformation.get_ligand_crds()));
		dbgmsg("ENERGY OF MINIMIZED LIGAND (DockedConformation) = " << min_conformation.get_energy());
		for (auto &pstate : min_conformation.get_states()) {
			dbgmsg("seed " << pstate->get_segment().get_seed_id() << " original E = " << pstate->get_energy() 
				<< " calculated E = " << __score.non_bonded_energy(__gridrec, pstate->get_segment().get_atoms(), pstate->get_crds()));
		}

#endif
                log_benchmark << "A-STAR for " << __ligand.name() << " took " << bench.seconds_from_start() 
                              << " wallclock seconds" << "\n";

                return __reconstruct(min_conformation);
        }

	Partial::Vec Linker::StaticLinker::__generate_rigid_conformations(const Seed::Graph &seed_graph) {
		Benchmark bench;
		log_note << "Generating rigid conformations of states for " << __ligand.name() << "..." << endl;

		State::Vec states;
		State::Id id = 0;
		for (auto &seed : seed_graph) {
#ifndef NDEBUG
			int state_no = 1;
#endif
			for (auto &pstate : seed.get_segment().get_states()) {
				states.push_back(&*pstate);
				pstate->set_id(id++);
#ifndef NDEBUG
				pstate->set_no(state_no++);
#endif
				dbgmsg(*pstate);
			}
		}
				
		//help::memusage("before find compatible state pairs");
	
		// init adjacency matrix (1 when states are compatible, 0 when not)
		Array2d<bool> conn = __find_compatible_state_pairs(seed_graph, states.size());

		dbgmsg("dimensions of conn = " << conn.get_szi() << " " << conn.get_szj());
		dbgmsg("conn = " << conn);
		//help::memusage("before max.clique.search");
		
		log_note << "find_compatible_state_pairs took " << bench.seconds_from_start() 
			<< " wallclock seconds for " << __ligand.name() << endl;
		
		Partial::Vec possibles_w_energy;
		{
			// find all maximum cliques, of size k; k...number of seed segments
			Maxclique m(conn);

			dbgmsg("largest possible clique (seed_graph.size()) = " << seed_graph.size());
			dbgmsg("currently searched for clique (__max_clique_size) = " << __max_clique_size);
			
			const int mcq_size = std::min(__max_clique_size,static_cast<int>(seed_graph.size()));
			const vector<vector<unsigned short int>> &qmaxes = m.mcq(mcq_size);
	
			//help::memusage("after max.clique.search");
	
			log_benchmark << "found " << qmaxes.size() << " maximum cliques, which took " 
				<< bench.seconds_from_start() << " wallclock seconds for " << __ligand.name() << "\n";
	
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
#ifndef NDEBUG
				stringstream ss;
				ss << "clique is composed of ";
				map<int, int> compo;
				
				for (auto &i : qmax) {
					ss << "seed_" << states[i]->get_segment().get_seed_id() 
						<< " state_" << states[i]->get_no() 
						<< " ene_" << states[i]->get_energy() << " ";
					compo[states[i]->get_segment().get_seed_id()] = states[i]->get_no();
				}
				ss << "ordered by seed_id : ";
				for (auto &kv : compo) {
					ss << kv.first << ":" << kv.second << " ";
				}
				ss << "total ene = " << energy;
				dbgmsg(ss.str());
#endif
			}
	
			//help::memusage("after possibles_w_energy");
		}
		//help::memusage("after forced destructor of qmaxes");
		
		// sort conformations according to their total energy
		Partial::sort(possibles_w_energy);

		if (static_cast<int>(possibles_w_energy.size()) > __max_num_possibles)
			possibles_w_energy.resize(__max_num_possibles);

		//help::memusage("after sort of possibles_w_energy");
		
		dbgmsg("number of possibles_w_energy = " << possibles_w_energy.size());
		// cluster rigid conformations and take only cluster representatives for further linking

		if (__docked_clus_rad <= 0.0) {
			log_benchmark << "Generated " << possibles_w_energy.size() 
				<< " possible conformations for ligand without clustering" << __ligand.name()
				<< ", which took " << bench.seconds_from_start() 
				<< " wallclock seconds" << "\n";
			return possibles_w_energy;
		}

		Partial::Vec clustered_possibles_w_energy = Cluster::Cluster::greedy(
			possibles_w_energy, __gridrec, __docked_clus_rad);
			
		//help::memusage("after greedy");
			
		dbgmsg("number of clustered_possibles_w_energy = " << clustered_possibles_w_energy.size());
		
		if (__max_possible_conf != -1 && static_cast<int>(clustered_possibles_w_energy.size()) > __max_possible_conf) {
			clustered_possibles_w_energy.resize(__max_possible_conf);
			dbgmsg("number of possible conformations > max_possible_conf, "
				<< "resizing to= " << __max_possible_conf << " conformations");
		}

		dbgmsg("RIGID CONFORMATIONS FOR LIGAND " << __ligand.name() 
			<< " : " << endl << clustered_possibles_w_energy);

		log_benchmark << "Generated " << clustered_possibles_w_energy.size() 
			<< " possible conformations for ligand " << __ligand.name()
			<< ", which took " << bench.seconds_from_start() 
			<< " wallclock seconds" << "\n";
		return clustered_possibles_w_energy;
	}
	
	DockedConformation Linker::StaticLinker::__reconstruct(const Partial &conformation) {
		Benchmark bench;
		log_note << "Reconstructing docked ligands for " << __ligand.name() << "..." << "\n";
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
			for (size_t i = 0; i < state.get_segment().get_atoms().size(); ++i) {
				Molib::Atom &atom = const_cast<Molib::Atom&>(state.get_segment().get_atom(i)); // ugly, correct this
				if (!overlap.count(&atom)) {
					const geometry::Coordinate &crd = state.get_crd(i);
					atom.set_crd(crd);
				}
			}
		}
		
		Molib::Molecule ligand(__ligand);
		ligand.set_name(__ligand.name() + "_" + std::to_string(++conf_number));
		
		log_benchmark << "Reconstruction of molecules took " << bench.seconds_from_start() 
			<< " wallclock seconds for " << __ligand.name() << "\n";
		
		// The receptor conformation will not be used, so don't copy it
		if (__max_iterations_final == -1)
			return DockedConformation(ligand, __receptor.name(), conformation.get_energy(), 0, 0);
		
		return DockedConformation(ligand, __receptor, conformation.get_energy(), 0, 0);
	}
};
}
