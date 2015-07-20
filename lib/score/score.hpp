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
	typedef map<const Atom*, Geom3D::Coordinate> AtomToCrd;
	class Score {
		typedef pair<int, int> pair_of_ints;
		typedef vector<double> M0;
		typedef map<pair_of_ints, M0> M1;
		M1 __gij_of_r_numerator;
		M1 __energies, __derivatives;
		map<pair_of_ints, double> __sum_gij_of_r_numerator;
		M0 __gij_of_r_bin_range_sum, __bin_range_sum;
		double __total_quantity;
		set<pair_of_ints> __prot_lig_pairs;
		const double __eps;
		const string __ref_state, __comp, __distributions_file, __rad_or_raw;
		const double __dist_cutoff, __step_non_bond;
		Atom::Grid &__gridrec;
		void __define_composition(const set<int>&, const set<int>&);
		void __process_distributions_file();
		void __compile_scoring_function();
		double __energy_mean(const pair_of_ints&, const double&);
		double __energy_cumulative(const pair_of_ints&, const double&);
		struct Inter {
			M0 potential, derivative;
		};
		Inter __interpolate(const M0&);
		int __get_index(const double d) const { return (int) floor(d / __step_non_bond); }
		double __get_lower_bound(const int idx) const { return (double) idx * __step_non_bond; }
	public:
		Score(const set<int> &receptor_idatm_types, const set<int> &ligand_idatm_types, 
				Atom::Grid &gridrec, const string &ref_state, 
				const string &comp, const string &rad_or_raw, const double &dist_cutoff, 
				const string &distributions_file, const double &step_non_bond) 
				: __gridrec(gridrec), __ref_state(ref_state), __comp(comp), 
				__rad_or_raw(rad_or_raw), __dist_cutoff(dist_cutoff), 
				__distributions_file(distributions_file), __step_non_bond(step_non_bond),
				__total_quantity(0), __eps(0.0000001) {
					
			function<double (Score&, const pair_of_ints&, const double&)> fptr = &Score::__energy_mean;
			__define_composition(receptor_idatm_types, ligand_idatm_types);
			__process_distributions_file();
			__compile_scoring_function();
		};
		double non_bonded_energy(const Molecule&) const; // this was formerly called distances_and_scores_frag_lig
		double non_bonded_energy(const Atom::Vec &atoms, const Geom3D::Point::Vec &crds) const;
		cluster::MapD<Molib::Molecule> many_ligands_score(const Molib::Molecules &ligands) const;
		Array1d<double> compute_energy(const Geom3D::Coordinate &crd, const set<int> &ligand_atom_types) const;
		const M1& get_energies() const { return __energies; }
		const M1& get_derivatives() const { return __derivatives; }
		friend ostream& operator<< (ostream& stream, const Score::M0 &energy);
	};
};
#endif
