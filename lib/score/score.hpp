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
		
		class InterpolationError : public Error {
		public: 
			InterpolationError(const string &msg) : Error(msg) {}
		};

		AtomPairValues __gij_of_r_numerator;

		AtomPairValues __energies, __derivatives; // objective function
		AtomPairValues __energies_scoring; // scoring function

		map<pair_of_ints, double> __sum_gij_of_r_numerator;
		vector<double> __gij_of_r_bin_range_sum, __bin_range_sum;
		double __total_quantity;
		set<pair_of_ints> __prot_lig_pairs;
		
		const double __eps;
		const string __ref_state, __comp, __rad_or_raw;
		const double __dist_cutoff, __step_non_bond;
		
		double __step_in_file;
		
		double __energy_mean(const pair_of_ints&, const double&);
		double __energy_cumulative(const pair_of_ints&, const double&);
		
		int __get_index(const double d) const { return (int) floor((d + 0.0000000001) / (double) __step_in_file); }
		double __get_lower_bound(const int idx) const { return (double) idx * (double) __step_in_file; }
	public:
		Score(const string &ref_state, const string &comp, const string &rad_or_raw, 
			const double &dist_cutoff, const double &step_non_bond) 
			: __total_quantity(0), __eps(0.0000001), __ref_state(ref_state), __comp(comp), __rad_or_raw(rad_or_raw),
			  __dist_cutoff(dist_cutoff), __step_non_bond(step_non_bond), __step_in_file(-1) {}
				
		double non_bonded_energy(const Atom::Grid &gridrec, const Molecule&) const; // this was formerly called distances_and_scores_frag_lig
		double non_bonded_energy(const Atom::Grid &gridrec, const Atom::Vec &atoms, const Geom3D::Point::Vec &crds) const;

		Array1d<double> compute_energy(const Atom::Grid &gridrec, const Geom3D::Coordinate &crd, const set<int> &ligand_atom_types) const;

		const AtomPairValues& get_energies() const { return __energies; }
		const AtomPairValues& get_derivatives() const { return __derivatives; }

		Score& define_composition(const set<int> &receptor_idatm_types, const set<int> &ligand_idatm_types);
		Score& process_distributions_file(const string &distributions_file);
		Score& compile_scoring_function();
		Score& compile_objective_function();
		Score& parse_objective_function(const string &obj_dir, const double scale_non_bond);
		Score& output_objective_function(const string &obj_dir);

		friend ostream& operator<< (ostream& stream, const vector<double> &energy);
		friend ostream& operator<< (ostream& stream, const Score::AtomPairValues &energies);
		friend ostream& operator<< (ostream& stream, const Score &score);
	};
};
#endif
