#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "molecule.hpp"
#include "geom3d/geom3d.hpp"
#include "fragmenter/fragmenter.hpp"
#include "fragmenter/unique.hpp"
#include "helper/benchmark.hpp"
#include "bond.hpp"
#include "bondtype.hpp"
using namespace std;

namespace Molib {
	ostream& operator<< (ostream& os, const ValenceState& valence_state) {
		for (auto &kv : valence_state) {
			const Atom &atom = *kv.first;
			const AtomParams &param = kv.second;
			os << "atom = " << atom.idatm_type_unmask() << "_" << atom.atom_number() 
				<< " valence = " << param.val << " aps = " << param.aps
				<< " connectivity = " << param.con << endl;
		}
		return os;
	}
	ostream& operator<< (ostream& os, const ValenceStateVec& valence_states) {
		for (int i = 0; i < valence_states.size(); ++i) {
			os << "VALENCE STATE NR. " << i << " : " << endl << valence_states[i];
		}
		return os;
	}
	ostream& operator<< (ostream& os, const BondToOrder& bond_orders) {
		for (auto &kv : bond_orders) {
			const Bond &b = *kv.first;
			//~ os << "bond_order(" << b.first_atom().atom_number() << ","
				//~ << b.second_atom().atom_number() << ") = " 
			os << "bond_order(" << b.atom1().atom_number() << ","
				<< b.atom2().atom_number() << ") = " 
				<< kv.second << endl;
		}
		return os;
	}
	//~ void BondOrder::compute_bond_order(Molecule &molecule) {
		//~ // assign atom types for atomic penalty scores
		//~ Fragmenter(molecule.get_atoms())
			//~ .substitute_atoms(help::atomic_penalty_scores);
		//~ dbgmsg("MOLECULE AFTER ASSIGNING ATOMIC PENALTY SCORES " << endl
			//~ << molecule);
		//~ // for tps 0, 1, 2, 3, ..., N create all possible combinations of valence states
		//~ // total number of valence states < 2000
		//~ const int max_valence_states = 2000;
		//~ ValenceStateVec valence_states = 
			//~ __create_valence_states(molecule, max_valence_states);
		//~ dbgmsg("VALENCE STATES ARE " << endl << valence_states);
		//~ // for each valence state determine bond orders (boaf procedure)
		//~ for (auto &valence_state : valence_states) {
			//~ try {
				//~ BondToOrder bond_orders;
				//~ if (__basic_rules(valence_state, bond_orders)) {
					//~ dbgmsg("successfully determined bond orders for molecule " 
						//~ << molecule.name() << " : " << endl << bond_orders);
					//~ // set the newly determined bond orders
					//~ for (auto &kv : bond_orders) {
						//~ Bond &bond = *kv.first;
						//~ bond.set_bo(kv.second);
					//~ }
					//~ return;
				//~ }
			//~ } catch(BondOrderError &e) {
				//~ cerr << e.what() << endl;
			//~ }
		//~ }
		//~ // if boaf fails for all saved valence states, a warning message is given
		//~ throw Error("die : bond order assignment failed for molecule "
			//~ + molecule.name());
	//~ }
	void BondOrder::compute_bond_gaff_type(Molecule &molecule) {
		// assign atom types for atomic penalty scores
		Fragmenter(molecule.get_atoms())
			.substitute_bonds(help::bond_gaff_type);
		dbgmsg("MOLECULE AFTER COMPUTING BOND GAFF TYPE " << endl
			<< molecule);
	}
	void BondOrder::compute_bond_order(Molecule &molecule) {
		try {
			// assign atom types for atomic penalty scores
			Fragmenter(molecule.get_atoms())
				.substitute_atoms(help::atomic_penalty_scores);
			dbgmsg("MOLECULE AFTER ASSIGNING ATOMIC PENALTY SCORES " << endl
				<< molecule);
			// for tps 0, 1, 2, 3, ..., N create all possible combinations of valence states
			// total number of valence states < 2000
			const int max_valence_states = 2000;
			ValenceStateVec valence_states = 
				__create_valence_states(molecule, max_valence_states);
			dbgmsg("VALENCE STATES ARE " << endl << valence_states << "-----");
			// for each valence state determine bond orders (boaf procedure)
			for (auto &valence_state : valence_states) {
				dbgmsg("determining bond orders for valence state");
				try {
					BondToOrder bond_orders;
					if (__basic_rules(valence_state, bond_orders)) {
						dbgmsg("successfully determined bond orders for molecule " 
							<< molecule.name() << " : " << endl << bond_orders);
						// set the newly determined bond orders
						for (auto &kv : bond_orders) {
							Bond &bond = *kv.first;
							bond.set_bo(kv.second);
						}
						dbgmsg("MOLECULE AFTER COMPUTING BOND ORDERS" 
							<< endl << molecule);
						return;
					}
				} catch(BondOrderError &e) {
					cerr << e.what() << endl;
				}
			}
			// if boaf fails for all saved valence states, a warning message is given
			throw Error("die : bond order assignment failed for molecule "
				+ molecule.name());
		} catch (exception& e) {
			cerr << "errmesg : " << e.what() << " for molecule = " << molecule.name() << endl;
		}
	}
	
	//~ void dfs(int level, int sum, int tps, const vector<vector<AtomParams> &V) {
		//~ auto &params = V[level];
		//~ for (auto &p : params) {
			//~ if (p.aps + sum < tps) {
				//~ Q.push(p);
				//~ dfs(level + 1, p.aps + sum, V[level + 1]);
				//~ Q.pop();
			//~ } else if (p.aps + sum == tps) {
				//~ // extend this valence state till the end
				//~ // save valence state
			//~ }
		//~ }
	//~ }
	void BondOrder::__dfs(const int level, const int sum, const int tps, 
		const vector<vector<AtomParams>> &V, vector<AtomParams> &Q, 
		vector<vector<AtomParams>> &valence_states, const int max_valence_states) {
		if (valence_states.size() > max_valence_states) return;
		auto &params = V[level];
		for (auto &p : params) {
			if (p.aps + sum <= tps) {
				Q.push_back(p);
				if (level + 1 < V.size()) {
					__dfs(level + 1, p.aps + sum, tps, V, Q,
						valence_states, max_valence_states);
				} else if (p.aps + sum == tps) {
					// save valence state from Q
					valence_states.push_back(Q);
				}
				Q.pop_back();
			}
		}
	}
	ValenceStateVec BondOrder::__create_valence_states(const Molecule &molecule, 
		const int max_valence_states) {
#ifndef NDEBUG
		Benchmark::reset();
#endif
		vector<Atom*> atoms;
		vector<vector<AtomParams>> V; // vector of AtomParams are sorted by increasing aps
		vector<vector<AtomParams>> valence_states;
		vector<AtomParams> Q;
		// prepare for recursive valence state determination
		for (auto &patom : molecule.get_atoms()) {
			Atom &atom = *patom;
			dbgmsg(atom);
			atoms.push_back(&atom);
			V.push_back(vector<AtomParams>());
			auto &aps = V.back();
			for (auto &kv : atom.get_aps()) {
				// sort increasingly atomic penalty scores for each atom
				aps.push_back(AtomParams{kv.first, kv.second, (int) atom.size()});
			}
			sort(aps.begin(), aps.end(), [] (const AtomParams &i, 
				const AtomParams &j) { return i.aps < j.aps; });
#ifndef NDEBUG			
			for (auto &x : aps) dbgmsg(x.val << " " << x.aps << " " << x.con);
#endif
		}
		// recursively find valence states
		for (int tps = 0; tps < 32; ++tps) {
			__dfs(0, 0, tps, V, Q, valence_states, max_valence_states);
		}
		// convert the result to atom mapping
		ValenceStateVec vss;
		for (auto &valence_state : valence_states) {
			vss.push_back(ValenceState());
			ValenceState &vs = vss.back();
			for (int i = 0; i < valence_state.size(); ++i) {
				vs[atoms[i]] = valence_state[i];
			}
		}
		dbgmsg(vss);
#ifndef NDEBUG
		double wall_secs = Benchmark::seconds_from_start();
		dbgmsg("Creating valence states took " << wall_secs << " seconds");
		Benchmark::reset();
#endif
		return vss;
	}

	//~ ValenceStateVec BondOrder::__create_valence_states(const Molecule &molecule, const int max_valence_states) {
		//~ ValenceStateVec valence_states;
		//~ AtomValencePairVec numbers;
		//~ for (auto &patom : molecule.get_atoms()) {
			//~ Atom &atom = *patom;
			//~ for (auto &kv : atom.get_aps()) {
				//~ numbers.push_back({&atom, 
					//~ AtomParams{kv.first, kv.second, (int) atom.size()}});
			//~ }
		//~ }
		//~ for (int tps = 0; tps < 1000; ++tps) {
			//~ if (!__sum_up_recursive(numbers, tps, AtomValencePairVec(), 
				//~ valence_states, max_valence_states))
				//~ break;
		//~ }
		//~ return valence_states;
	//~ }
	//~ bool BondOrder::__sum_up_recursive(const AtomValencePairVec &numbers, const int tps, 
		//~ const AtomValencePairVec &partial, ValenceStateVec &valence_states,
		//~ const int max_valence_states) {
		//~ if (valence_states.size() <= max_valence_states) {
			//~ int s = 0;
			//~ for (auto &x : partial) 
				//~ s += x.second.aps;
			//~ if (s == tps) {
				//~ ValenceState valence_state;
				//~ for (auto &pp : partial)
					//~ valence_state.insert({pp.first, pp.second});
				//~ valence_states.push_back(valence_state);
			//~ } else if (s > tps)
				//~ return true;
			//~ for (int i = 0; i < numbers.size(); ++i) {
				//~ AtomValencePairVec remaining;
				//~ auto n = numbers[i];
				//~ for (int j = i + 1; j < numbers.size(); ++j) {
					//~ if (numbers[i].first != numbers[j].first) {
						//~ remaining.push_back(numbers[j]);
					//~ }
				//~ }
				//~ AtomValencePairVec partial_rec = partial;
				//~ partial_rec.push_back(n);
				//~ if (!__sum_up_recursive(remaining, tps, partial_rec, 
					//~ valence_states, max_valence_states))
					//~ return false;
			//~ }
		//~ }
		//~ dbgmsg("max number of valence states (" << max_valence_states 
			//~ << ") exceeded at tps = " << tps);
		//~ return false;
	//~ }
	bool BondOrder::__discrepancy(const ValenceState &valence_state) {
		// b) if discrepancy happens (av is not 0 when con is 0 or av is 0 when con is not 0)
		// reset the bond order to 2 and then 3.
		for (auto &kv : valence_state) {
			const Atom &atom = *kv.first;
			const AtomParams &apar = kv.second;
			if (apar.val != 0 && apar.con == 0
				|| apar.val == 0 && apar.con != 0) {
				return true;
			}
		}
		return false;
	}
	bool BondOrder::__success(const ValenceState &valence_state) {
		// rule 4 : if all the bonds are successfully assigned, con and av 
		// of every atom are both 0 (boaf returns 1 and stops)
		for (auto &kv : valence_state) {
			const AtomParams &apar = kv.second;
			if (apar.val != 0 || apar.con != 0) {
				return false;
			}
		}
		return true;
	}
	bool BondOrder::__basic_rules(ValenceState &valence_state, BondToOrder &bond_orders) {
		while (!__success(valence_state)) {
			bool bo_was_set = false;
			for (auto &kv : valence_state) {
				Atom &atom = *kv.first;
				AtomParams &apar = kv.second;
				// rule 2 : for one atom, if its con equals to av, the bond 
				// orders of its unassigned bonds are set to 1
				if (apar.con == apar.val) {
					//~ for (auto &bond : atom)
					for (auto &pbond : atom.get_bonds()) {
						Bond &bond = *pbond;
						if (!bond_orders.count(&bond)) { // unassigned bond order
							// rule 1 : for each atom in a bond, if the bond order bo 
							// is determined, con is deducted by 1 and av is deducted by bo
							bond_orders[&bond] = 1;
							//~ bond_orders[&bond.get_reverse()] = 1;
							//~ AtomParams &apar2 = valence_state.at(&bond.second_atom());
							AtomParams &apar2 = valence_state.at(&bond.second_atom(atom));
							apar.con -= 1;
							apar.val -= 1;
							apar2.con -= 1;
							apar2.val -= 1;
							bo_was_set = true;
						}
					}
				}
				// rule 3 : for one atom, if its con equals to 1, the bond 
				// order of the last bond is set to av
				else if (apar.con == 1) {
					//~ for (auto &bond : atom)
					for (auto &pbond : atom.get_bonds()) {
						Bond &bond = *pbond;
						if (!bond_orders.count(&bond)) { // unassigned bond order
							// rule 1 : for each atom in a bond, if the bond order bo 
							// is determined, con is deducted by 1 and av is deducted by bo
							const int bo = apar.val;
							bond_orders[&bond] = bo;
							//~ bond_orders[&bond.get_reverse()] = bo;
							//~ AtomParams &apar2 = valence_state.at(&bond.second_atom());
							AtomParams &apar2 = valence_state.at(&bond.second_atom(atom));
							apar.con -= 1;
							apar.val -= bo;
							apar2.con -= 1;
							apar2.val -= bo;
							bo_was_set = true;
						}
					}
				}
			}
			// b) if discrepancy happens (av is not 0 when con is 0 or av is 0 when con is not 0)
			// reset the bond order to 2 and then 3.
			if (__discrepancy(valence_state))
				return false;
			// if non of the above rules can be applied do the trial-error test :
			if (!bo_was_set && !__success(valence_state)) 
				__trial_error(valence_state, bond_orders);
		}
		return true;
	}
	Bond& BondOrder::__get_first_unassigned_bond(const ValenceState &valence_state, BondToOrder &bond_orders) {
		for (auto &kv : valence_state) {
			const Atom &atom = *kv.first;
			const AtomParams &apar = kv.second;
			//~ for (auto &bond : atom)
			for (auto &pbond : atom.get_bonds()) {
				Bond &bond = *pbond;
				if (!bond_orders.count(&bond)) { // unassigned bond order
					return bond;
				}
			}
		}
		throw BondOrderError("exception : cannot find unassigned bond");
	}

	void BondOrder::__trial_error(ValenceState &valence_state, BondToOrder &bond_orders) {
		for (int bo = 1; bo <= 3; ++bo) {
			// save valence state
			ValenceState saved = valence_state;
			BondToOrder saved_bond_orders = bond_orders;
			// a) assume bond order (of one randomly selected bond) is 1, then continue 
			// bond order assignemnt according to the basic rules
			Bond &bond = __get_first_unassigned_bond(valence_state, bond_orders);
			bond_orders[&bond] = bo;
			//~ bond_orders[&bond.get_reverse()] = bo;
			//~ AtomParams &apar1 = valence_state.at(&bond.first_atom());
			//~ AtomParams &apar2 = valence_state.at(&bond.second_atom());
			AtomParams &apar1 = valence_state.at(&bond.atom1());
			AtomParams &apar2 = valence_state.at(&bond.atom2());
			apar1.con -= 1;
			apar1.val -= bo;
			apar2.con -= 1;
			apar2.val -= bo;
			if (__basic_rules(valence_state, bond_orders)) return;
			// reset valence state to the saved one
			valence_state = saved;
			bond_orders = saved_bond_orders;
		}
		// c) if discrepancies happen for all 3 bond orders, boaf exits and returns 0	
		throw BondOrderError("exception : discrepancies happened for all 3 bond orders");	
	}
};
