#ifndef SCORE_H
#define SCORE_H
#include "candock/helper/benchmark.hpp"
#include "candock/helper/error.hpp"
#include "candock/molib/atom.hpp"
#include "candock/cluster/optics.hpp"
#include "candock/geometry/geometry.hpp"
#include <iostream>
#include <exception>
#include <typeinfo>
#include <map>
#include <set>
#include <cmath>
#include "candock/helper/array1d.hpp"

namespace candock {

namespace molib {
        class Molecule;
        class Molecules;
}

namespace score {

        class Score {
        protected:
                typedef std::pair<int, int> pair_of_ints;
                typedef std::map<pair_of_ints, std::vector<double>> AtomPairValues;

                AtomPairValues __gij_of_r_numerator;

                AtomPairValues __energies_scoring; // scoring function

                std::map<pair_of_ints, double> __sum_gij_of_r_numerator;
                std::vector<double> __gij_of_r_bin_range_sum, __bin_range_sum;
                double __total_quantity;
                std::set<pair_of_ints> __prot_lig_pairs;
                std::set<pair_of_ints> __avail_prot_lig;

                const double __eps;
                const std::string __ref_state, __comp, __rad_or_raw;
                const double __dist_cutoff;

                double __step_in_file;

                double __energy_mean(const pair_of_ints&, const double&);
                double __energy_cumulative(const pair_of_ints&, const double&);

                int __get_index(const double d) const { return (int) floor((d + 0.0000000001) / (double) __step_in_file); }
                double __get_lower_bound(const int idx) const { return (double) idx * (double) __step_in_file; }
        public:
                Score(const std::string &ref_state, const std::string &comp, const std::string &rad_or_raw, 
                        const double &dist_cutoff) 
                        : __total_quantity(0), __eps(0.0000001), __ref_state(ref_state), __comp(comp), __rad_or_raw(rad_or_raw),
                          __dist_cutoff(dist_cutoff), __step_in_file(-1) {}

                double non_bonded_energy(const molib::Atom::Grid &gridrec, const molib::Molecule&) const; // this was formerly called distances_and_scores_frag_lig
                double non_bonded_energy(const molib::Atom::Grid &gridrec, const molib::Atom::Vec &atoms, const geometry::Point::Vec &crds) const;

                Array1d<double> compute_energy(const molib::Atom::Grid &gridrec, const geometry::Coordinate &crd, const std::set<int> &ligand_atom_types) const;

                double get_dist_cutoff() const { return __dist_cutoff; }

                Score& define_composition(const std::set<int> &receptor_idatm_types, const std::set<int> &ligand_idatm_types);
                Score& process_distributions_file(const std::string &distributions_file);
                Score& compile_scoring_function();

                friend std::ostream& operator<< (std::ostream& stream, const std::vector<double> &energy);
                friend std::ostream& operator<< (std::ostream& stream, const Score::AtomPairValues &energies);
                friend std::ostream& operator<< (std::ostream& stream, const Score &score);
        };
}

}

#endif