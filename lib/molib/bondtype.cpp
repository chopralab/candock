#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <numeric>  
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "candock/molib/molecule.hpp"
#include "candock/geom3d/geom3d.hpp"
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/fragmenter/unique.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/renamerules.hpp"
#include "candock/molib/bond.hpp"
#include "candock/molib/bondtype.hpp"
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
		for (size_t i = 0; i < valence_states.size(); ++i) {
			os << "VALENCE STATE NR. " << i << " : " << endl << valence_states[i];
		}
		return os;
	}
	ostream& operator<< (ostream& os, const BondToOrder& bond_orders) {
		for (auto &kv : bond_orders) {
			const Bond &b = *kv.first;
			os << "bond_order(" << b.atom1().atom_number() << ","
				<< b.atom2().atom_number() << ") = " 
				<< kv.second << endl;
		}
		return os;
	}
	
	void BondOrder::compute_rotatable_bonds(const Atom::Vec &atoms) {
		dbgmsg("starting compute_rotatable_bonds");
		Fragmenter f(atoms);
		f.substitute_bonds(help::rotatable);
		// rotatable bond inside a ring is changed back to non-rotatable
		Rings rings = f.identify_fused_rings();
		//~ Rings rings = f.identify_rings();
		for (auto &ring : rings) {
			for (auto &pbond : get_bonds_in(ring)) {
				pbond->set_rotatable("");
			}
		}
		dbgmsg("MOLECULE AFTER COMPUTING ROTATABLE BONDS" << endl << atoms);
	}

	
	void BondOrder::compute_bond_gaff_type(const Atom::Vec &atoms) {
		// assign atom types for atomic penalty scores
		Fragmenter(atoms)
			.substitute_bonds(help::bond_gaff_type);
		dbgmsg("MOLECULE AFTER COMPUTING BOND GAFF TYPE " << endl << atoms);
	}
	void BondOrder::compute_bond_order(const Atom::Vec &atoms) {
		try {
			// assign atom types for atomic penalty scores
			Fragmenter(atoms)
				.substitute_atoms(help::atomic_penalty_scores);
			dbgmsg("MOLECULE AFTER ASSIGNING ATOMIC PENALTY SCORES " << endl << atoms);
			// for tps 0, 1, 2, 3, ..., N create all possible combinations of valence states
			// total number of valence states < 2000
			const int max_valence_states = 2000;
			ValenceStateVec valence_states = 
				__create_valence_states(atoms, max_valence_states);
			dbgmsg("VALENCE STATES ARE " << endl << valence_states << "-----");
			// for each valence state determine bond orders (boaf procedure)
#ifndef NDEBUG
			int val_cnt = 0;
#endif
			for (auto &valence_state : valence_states) {
				dbgmsg("determining bond orders for valence state " << val_cnt++);
				try {
					BondToOrder bond_orders;
					if (__basic_rules(valence_state, bond_orders)) {
						dbgmsg("successfully determined bond orders for molecule : " << endl << bond_orders);
						// set the newly determined bond orders
						for (auto &kv : bond_orders) {
							Bond &bond = *kv.first;
							bond.set_bo(kv.second);
						}
						dbgmsg("MOLECULE AFTER COMPUTING BOND ORDERS" << endl << atoms);
						return;
					}
				} catch(BondOrderError &e) {
					dbgmsg(e.what());
				}
			}
			// if boaf fails for all saved valence states, a warning message is given
			throw Error("die : bond order assignment failed");
		} catch (exception& e) {
			log_error << "errmesg : " << e.what() << " for residue = " << endl << atoms << endl;
			throw e;
		}
	}

        void compute_chirality(const Atom::Vec &bonds) {
                for (auto a : bonds) {
                        if (a->idatm_type_unmask() != "C3") {
                                continue;
                        }
                        if (a->get_num_hydrogens() != 1) {
                                continue;
                        }

                        auto atom_bonds = a->get_bonds();

                        if (atom_bonds.size() != 4) {
                            throw BondOrderError("Error assigning chirality, wrong number of neighbors");
                        }

                        // FIXME: Using the topologic order is a bad idea, but it's all I got for now
                        Atom::Vec neighbors;
                        for (auto b : atom_bonds) {
                                neighbors.push_back(&b->second_atom(*a));
                        }

                        atom_bonds[0]->set_stereo("I");
                        atom_bonds[1]->set_stereo("I");

                        double improper_1 = Geom3D::dihedral(neighbors[0]->crd(), a->crd(), neighbors[1]->crd(), neighbors[2]->crd());
                        double improper_2 = Geom3D::dihedral(neighbors[0]->crd(), a->crd(), neighbors[1]->crd(), neighbors[3]->crd());

                        if (improper_1 > improper_2) {
                                atom_bonds[2]->set_stereo("U");
                                atom_bonds[3]->set_stereo("D");
                        } else {
                                atom_bonds[2]->set_stereo("D");
                                atom_bonds[3]->set_stereo("U");
                        }
                }
        }

	void BondOrder::__dfs(const int level, const int sum, const int tps, 
		const vector<vector<AtomParams>> &V, vector<AtomParams> &Q, 
		vector<vector<AtomParams>> &valence_states, const int max_valence_states) {
		if (static_cast<int>(valence_states.size()) > max_valence_states) return;
		auto &params = V[level];
		for (auto &p : params) {
			if (p.aps + sum <= tps) {
				Q.push_back(p);
				if (level + 1 < static_cast<int>(V.size())) {
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

	ValenceStateVec BondOrder::__create_valence_states(const Atom::Vec &atoms, const int max_valence_states) {
#ifndef NDEBUG
		Benchmark bench;
#endif
		vector<vector<AtomParams>> V; // vector of AtomParams are sorted by increasing aps
		vector<vector<AtomParams>> valence_states;
		vector<AtomParams> Q;
		// prepare for recursive valence state determination
		for (auto &patom : atoms) {
			Atom &atom = *patom;
			dbgmsg(atom);
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

		// sort valence states within each tps from lower to higher individual aps (issue #103)
		sort(valence_states.begin(), valence_states.end(), [] (const vector<AtomParams> &i, const vector<AtomParams> &j) { 
				int tps_i = 0; for (auto &e : i) tps_i += e.aps;
				int tps_j = 0; for (auto &e : j) tps_j += e.aps;
				const AtomParams max_ap_i = *max_element(i.begin(), i.end(), [] (const AtomParams &iap, const AtomParams &jap) { return iap.aps < jap.aps; });
				const AtomParams max_ap_j = *max_element(j.begin(), j.end(), [] (const AtomParams &iap, const AtomParams &jap) { return iap.aps < jap.aps; });
				if (tps_i == tps_j)
					return max_ap_i.aps < max_ap_j.aps;
				else
					return tps_i < tps_j;
			});
		
		// convert the result to atom mapping
		ValenceStateVec vss;
		for (auto &valence_state : valence_states) {
			vss.push_back(ValenceState());
			ValenceState &vs = vss.back();
			for (size_t i = 0; i < valence_state.size(); ++i) {
				vs[atoms[i]] = valence_state[i];
			}
		}
		dbgmsg(vss);
#ifndef NDEBUG
		double wall_secs = bench.seconds_from_start();
		dbgmsg("Creating valence states took " << wall_secs << " seconds");
		bench.reset();
#endif
		return vss;
	}

	bool BondOrder::__discrepancy(const ValenceState &valence_state) {
		// b) if discrepancy happens (av is not 0 when con is 0 or av is 0 when con is not 0)
		// reset the bond order to 2 and then 3.
		for (auto &kv : valence_state) {
			const AtomParams &apar = kv.second;
			if ( (apar.val != 0 && apar.con == 0)
                          || (apar.val == 0 && apar.con != 0)) {
                            dbgmsg( "Failure for atom = " << kv.first->atom_name() << " val = " << apar.val
                                << " con = " << apar.con << " aps = " << apar.aps);
                                return true;
			}
		}
		return false;
	}

	bool BondOrder::__my_success(const ValenceState &valence_state) {
		// rule 4 : if all the bonds are successfully assigned, con and av 
		// of every atom are both 0 (boaf returns 1 and stops)
		for (auto &kv : valence_state) {
			const AtomParams &apar = kv.second;
			dbgmsg("atom = " << kv.first->atom_name() << " val = " << apar.val
				<< " con = " << apar.con << " aps = " << apar.aps);
			if (apar.val != 0 || apar.con != 0) {
				return false;
			}
		}
		return true;
	}

	bool BondOrder::__basic_rules(ValenceState &valence_state, BondToOrder &bond_orders) {
		while (!__my_success(valence_state)) {
			bool bo_was_set = false;
			for (auto &kv : valence_state) {
				Atom &atom = *kv.first;
				AtomParams &apar = kv.second;
				// rule 2 : for one atom, if its con equals to av, the bond 
				// orders of its unassigned bonds are set to 1
				if (apar.con == apar.val) {
					for (auto &pbond : atom.get_bonds()) {
						Bond &bond = *pbond;
						if (!bond_orders.count(&bond)) { // unassigned bond order
							// rule 1 : for each atom in a bond, if the bond order bo 
							// is determined, con is deducted by 1 and av is deducted by bo
							bond_orders[&bond] = 1;
                                                        if (valence_state.count(&bond.second_atom(atom)) == 0 ) {
                                                                apar.con -= 1;
                                                                apar.val -= 1;
                                                                bo_was_set = true;
                                                                continue;
                                                        }
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
					for (auto &pbond : atom.get_bonds()) {
						Bond &bond = *pbond;
						if (!bond_orders.count(&bond)) { // unassigned bond order
							// rule 1 : for each atom in a bond, if the bond order bo 
							// is determined, con is deducted by 1 and av is deducted by bo
							const int bo = apar.val;
							// rule NEW (Janez) : if trying to assign a zero bond order then 
							// try another valence state
							dbgmsg("bo =" << bo);
							if (bo <= 0)
								throw BondOrderError("exception : zero bond order");
							bond_orders[&bond] = bo;
                                                        if (valence_state.count(&bond.second_atom(atom)) == 0 ) {
                                                                apar.con -= 1;
                                                                apar.val -= bo;
                                                                bo_was_set = true;
                                                                continue;
                                                        }
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
			if (!bo_was_set && !__my_success(valence_state)) 
				__trial_error(valence_state, bond_orders);
		}
		return true;
	}
	Bond& BondOrder::__get_first_unassigned_bond(const ValenceState &valence_state, BondToOrder &bond_orders) {
		for (auto &kv : valence_state) {
			const Atom &atom = *kv.first;
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
