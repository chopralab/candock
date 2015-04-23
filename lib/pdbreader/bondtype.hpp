#ifndef BONDTYPE_H
#define BONDTYPE_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <cstdlib>
#include "it.hpp"

using namespace std;

namespace Molib {	
	class Atom;
	struct AtomParams {
		int val;
		int aps;
		int con;
	};
	typedef map<Atom*, AtomParams> ValenceState;
	typedef vector<ValenceState> ValenceStateVec;
	typedef map<Bond*, int> BondToOrder;
	
	class BondOrder {
		class BondOrderError : public Error {
		public: 
			BondOrderError(const string &msg) : Error(msg) {}
		};
		static ValenceStateVec __create_valence_states(const Molecule &molecule, const int max_valence_states);
		static void __dfs(const int level, const int sum, const int tps, const vector<vector<AtomParams>> &V,
			vector<AtomParams> &Q, vector<vector<AtomParams>> &valence_states, const int max_valence_states);
		static bool __discrepancy(const ValenceState &valence_state);
		static bool __success(const ValenceState &valence_state);
		static bool __basic_rules(ValenceState &valence_state, BondToOrder &bond_orders);
		static void __trial_error(ValenceState &valence_state, BondToOrder &bond_orders);
		static Bond& __get_first_unassigned_bond(const ValenceState &valence_state, BondToOrder &bond_orders);
	public:
		static void compute_bond_order(Molecule &molecule);
		static void compute_bond_gaff_type(Molecule &molecule);
	};
	ostream& operator<< (ostream& stream, const ValenceState& valence_state);
	ostream& operator<< (ostream& stream, const ValenceStateVec& valence_states);
	ostream& operator<< (ostream& os, const BondToOrder& bond_orders);
};
#endif
