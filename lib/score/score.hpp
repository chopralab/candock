#ifndef SCORE_H
#define SCORE_H
#include "helper/benchmark.hpp"
#include "helper/error.hpp"
#include "pdbreader/atom.hpp"
#include "cluster/optics.hpp"
#include "geom3d/geom3d.hpp"
#include <iostream>
#include <exception>
#include <typeinfo>
#include <map>
#include <set>
#include <cmath>
#include "helper/array1d.hpp"
using namespace std;

namespace Molib {

	class Molecule;
	class Molecules;

	class Score {
		typedef pair<int, int> pair_of_ints;
		typedef map<pair_of_ints, vector<double>> AtomPairValues;
		AtomPairValues __gij_of_r_numerator;

		AtomPairValues __energies, __derivatives; // objective function
		AtomPairValues __energies_scoring; // scoring function

		map<pair_of_ints, double> __sum_gij_of_r_numerator;
		vector<double> __gij_of_r_bin_range_sum, __bin_range_sum;
		double __total_quantity;
		set<pair_of_ints> __prot_lig_pairs;
		
		const double __eps;
		const string __ref_state, __comp, __distributions_file, __rad_or_raw;
		const double __dist_cutoff, __step_non_bond, __scale_non_bond;
		
		double __step_in_file;
		
		void __define_composition(const set<int>&, const set<int>&);
		void __process_distributions_file();
		void __compile_objective_function();
		void __compile_scoring_function();
		double __energy_mean(const pair_of_ints&, const double&);
		double __energy_cumulative(const pair_of_ints&, const double&);
		
		int __get_index(const double d) const { return (int) floor(d + 0.0001 / (const double) __step_in_file); }
		double __get_lower_bound(const int idx) const { return (double) idx * (const double) __step_in_file; }
	public:
		Score(const set<int> &receptor_idatm_types, const set<int> &ligand_idatm_types, const string &ref_state, 
				const string &comp, const string &rad_or_raw, const double &dist_cutoff, 
				const string &distributions_file, const double &step_non_bond, const double &scale_non_bond) 
				: __ref_state(ref_state), __comp(comp), 
				__rad_or_raw(rad_or_raw), __dist_cutoff(dist_cutoff), 
				__distributions_file(distributions_file), __step_non_bond(step_non_bond),
				__scale_non_bond(scale_non_bond), __total_quantity(0), __eps(0.0000001), __step_in_file(-1) {
			
			try {
				function<double (Score&, const pair_of_ints&, const double&)> fptr = &Score::__energy_mean;
				__define_composition(receptor_idatm_types, ligand_idatm_types);
				__process_distributions_file();
				__compile_objective_function();
				__compile_scoring_function();
			} catch (Error &e) {
				cerr << "Error in constructor of Score : " << e.what() << endl;
				throw e;
			}
		};
		double non_bonded_energy(Atom::Grid &gridrec, const Molecule&) const; // this was formerly called distances_and_scores_frag_lig
		double non_bonded_energy(Atom::Grid &gridrec, const Atom::Vec &atoms, const Geom3D::Point::Vec &crds) const;

		Array1d<double> compute_energy(Atom::Grid &gridrec, const Geom3D::Coordinate &crd, const set<int> &ligand_atom_types) const;

		friend ostream& operator<< (ostream& stream, const vector<double> &energy);
		friend ostream& operator<< (ostream& stream, const Score::AtomPairValues &energies);
		friend ostream& operator<< (ostream& stream, const Score &score);
	};
};
#endif
