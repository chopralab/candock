#ifndef KBFF_H
#define KBFF_H
#include "candock/score/score.hpp"

namespace Score {
    class KBFF : public Score {
        AtomPairValues __energies, __derivatives; // objective function
        double __step_non_bond;
    public:
        KBFF(const std::string &ref_state, const std::string &comp, const std::string &rad_or_raw, 
              const double &dist_cutoff, const double &step_non_bond) :
            Score(ref_state, comp, rad_or_raw, dist_cutoff), __step_non_bond(step_non_bond)
        {}

        double get_step_nonbond()const { return __step_non_bond; }
        const AtomPairValues& get_energies() const { return __energies; }
        const AtomPairValues& get_derivatives() const { return __derivatives; }

        KBFF& compile_objective_function(const double scale_non_bond);
        KBFF& parse_objective_function(const std::string &obj_dir, const double scale_non_bond, const size_t max_step);
        KBFF& output_objective_function(const std::string &obj_dir);

        friend ostream& operator<< (ostream& stream, const KBFF &kbff);
    };
}

#endif