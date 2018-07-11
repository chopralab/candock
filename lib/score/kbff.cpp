#include "candock/score/kbff.hpp"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include "candock/helper/help.hpp"
#include "candock/helper/path.hpp"
#include "candock/score/interpolation.hpp"
#include "candock/score/powerfit.hpp"
using namespace std;

namespace candock {

namespace score {

ostream& operator<<(ostream& stream, const KBFF& score) {
    for (auto& kv : score.__energies_scoring) {
        auto& atom_pair = kv.first;
        auto& ene = score.__energies.at(atom_pair);
        auto& der = score.__derivatives.at(atom_pair);
        stream << "#objective function" << endl;
        for (size_t i = 0; i < ene.size(); ++i) {
            stream << help::idatm_unmask[atom_pair.first] << "\t"
                   << help::idatm_unmask[atom_pair.second] << "\t" << fixed
                   << setprecision(3) << i * score.__step_non_bond << "\t"
                   << fixed << setprecision(3) << ene[i] << "\t" << fixed
                   << setprecision(3) << der[i] << endl;
        }
    }
    return stream;
}

KBFF& KBFF::output_objective_function(const string& obj_dir) {
    for (auto& kv : __energies) {
        auto& atom_pair = kv.first;
        auto& ene = __energies.at(atom_pair);
        //~ auto &der = __derivatives.at(atom_pair);

        stringstream ss;

        const string idatm_type1 = help::idatm_unmask[atom_pair.first];
        const string idatm_type2 = help::idatm_unmask[atom_pair.second];

        for (size_t i = 0; i < ene.size(); ++i) {
            ss << setprecision(8) << ene[i] << endl;
        }
        const string& filename = idatm_type1 + "_" + idatm_type2 + ".txt";

        Inout::file_open_put_stream(
            Path::join(Path::join(obj_dir, std::to_string(__step_non_bond)),
                       filename),
            ss);
    }
    return *this;
}

KBFF& KBFF::parse_objective_function(const string& obj_dir,
                                     const double scale_non_bond,
                                     const size_t max_step) {
    log_step << "Parsing obj function from disk" << endl;

    boost::filesystem::path path_to_objective_function(obj_dir);
    std::string subdir = std::to_string(__step_non_bond);

    // Some systems add trailing zeros to doubles, let's remove them
    while (!boost::filesystem::exists(path_to_objective_function / subdir) &&
           subdir.back() == '0') {
        subdir.erase(subdir.end() - 1);
    }

    path_to_objective_function /= subdir;

    if (!boost::filesystem::exists(path_to_objective_function)) {
        throw Error("Objective function not found! check 'obj_dir' and 'step'");
    }

    for (auto& atom_pair : __avail_prot_lig) {
        const string idatm_type1 = help::idatm_unmask[atom_pair.first];
        const string idatm_type2 = help::idatm_unmask[atom_pair.second];

        dbgmsg("parsing objective function for " << idatm_type1 << " and "
                                                 << idatm_type2);

        const string& filename = idatm_type1 + "_" + idatm_type2 + ".txt";

        // Check if the file exists on disk. Since we know that the directory
        // exists on disk,
        // we can skip this pair as it is not part of the scoring function (ie
        // not in the CSD).
        if (!Inout::file_size(
                (path_to_objective_function / filename).string())) {
            continue;
        }

        vector<string> contents;
        Inout::read_file((path_to_objective_function / filename).string(),
                         contents);

        for (auto& line : contents) {
            stringstream ss(line);
            string str1;
            ss >> str1;
            __energies[atom_pair].push_back(stod(str1));
        }

        if (__energies[atom_pair].size() < max_step) {
            int energies_diff = max_step - __energies[atom_pair].size();

            if (energies_diff < 0) {
                log_warning
                    << "Large atom pair: " << __energies[atom_pair].size()
                    << " " << idatm_type1 << " " << idatm_type2 << endl;
                __energies[atom_pair].resize(max_step);
            }

            for (int i = 0; i < energies_diff; ++i) {
                __energies[atom_pair].push_back(0.0);
            }
        }

        __derivatives[atom_pair] =
            Interpolation::derivative(__energies[atom_pair], __step_non_bond);

        // scale derivatives ONLY not energies and don't do it in kbforce cause
        // here is more efficient
        for (auto& dEdt : __derivatives[atom_pair]) {
            dEdt *= scale_non_bond;
        }
    }
    dbgmsg("parsed objective function");
    return *this;
}

KBFF& KBFF::compile_objective_function(const double scale_non_bond) {
    log_step << "Compiling objective function for minimization...\n";
    Benchmark bench;
    auto energy_function = __ref_state == "mean"
                               ? mem_fn(&KBFF::__energy_mean)
                               : mem_fn(&KBFF::__energy_cumulative);
    for (auto& el1 : __gij_of_r_numerator) {
        const pair_of_ints& atom_pair = el1.first;

        // Don't bother making an objective function
        // unless we're going to use it
        if (!__avail_prot_lig.count(atom_pair)) continue;

        const vector<double>& gij_of_r_vals = el1.second;
        dbgmsg(atom_pair.first << " " << atom_pair.second);
        const double w1 = help::vdw_radius[atom_pair.first];
        const double w2 = help::vdw_radius[atom_pair.second];
        const double vdW_sum = ((w1 > 0 && w2 > 0) ? w1 + w2 : 4.500);

        const string idatm_type1 = help::idatm_unmask[atom_pair.first];
        const string idatm_type2 = help::idatm_unmask[atom_pair.second];

        vector<double> energy(gij_of_r_vals.size(), -HUGE_VAL);

        const size_t start_idx = __get_index(vdW_sum - 0.6);
        const size_t end_idx = __get_index(vdW_sum + 1.0);

        dbgmsg(start_idx << " " << end_idx);

        for (size_t i = 0; i < gij_of_r_vals.size(); ++i) {
            const double lower_bound = __get_lower_bound(i);
            dbgmsg(lower_bound);
            const double& gij_of_r_numerator = gij_of_r_vals[i];
            dbgmsg("lower bound for atom_pair "
                   << idatm_type1 << " " << idatm_type2 << " " << lower_bound
                   << " gij_of_r_numerator = " << gij_of_r_numerator
                   << " __sum_gij_of_r_numerator[atom_pair] = "
                   << __sum_gij_of_r_numerator[atom_pair]);

            if (__sum_gij_of_r_numerator[atom_pair] < __eps) {
                energy[i] = 0;
            } else if (gij_of_r_numerator < __eps && (i + 1 < start_idx)) {
                energy[i] = 5.0;
            } else if (gij_of_r_numerator < __eps && (i + 1 >= start_idx)) {
                energy[i] = 0;
            } else {
                energy[i] = energy_function(*this, atom_pair, lower_bound);
            }
        }

        // dbgmsg("raw energies before interpolations : " << std::endl <<
        // energy);

        // correct for outliers only in the TRUE potential region of the
        // potential
        for (size_t i = start_idx; i < energy.size() - 2; ++i) {
            if (energy[i] != 0 && energy[i] != 5) {
                size_t j = i + 1;
                while (j < i + 3 && j < energy.size() &&
                       (energy[j] == 0 || energy[j] == 5)) {
                    ++j;
                }
                if (j > i + 1 && j < energy.size()) {
                    const double k0 = (energy[j] - energy[i]) / (j - i);
                    for (size_t k = 1; k < j - i; ++k) {
                        energy[i + k] = energy[i] + k0 * k;
                    }
                }
            }
        }

        // locate global minimum and repulsion index in interval [vdW_sum - 0.6,
        // vdW_sum + 1.0]
        double global_min = HUGE_VAL;
        size_t global_min_idx = start_idx;
        size_t repulsion_idx = global_min_idx;

        try {
            repulsion_idx = __get_index(
                help::repulsion_idx.at(make_pair(idatm_type1, idatm_type2)));
        } catch (out_of_range&) {
            try {
                repulsion_idx = __get_index(help::repulsion_idx.at(
                    make_pair(idatm_type2, idatm_type1)));
            } catch (out_of_range&) {
                dbgmsg("de-novo calculation of repulsion idx");

                for (size_t i = start_idx; i < end_idx && i < energy.size();
                     ++i) {
                    if (energy[i] != -HUGE_VAL) {
                        if (energy[i] < global_min) {
                            global_min = energy[i];
                            global_min_idx = i;

                            // minimum has to have steep downward slope on the
                            // left side
                            // leave the slope intact and unset everything left
                            // of the slope start point
                            repulsion_idx = i;
                            while (repulsion_idx > 0 &&
                                   energy[repulsion_idx] != -HUGE_VAL &&
                                   ((energy[repulsion_idx] -
                                     energy[repulsion_idx - 1]) /
                                    __step_in_file) < 0.75)
                                --repulsion_idx;
                        }
                    }
                }
            }
        }

        dbgmsg("repulsion idx (before correction) = "
               << repulsion_idx
               << " repulsion distance (below is forbidden area) = "
               << __get_lower_bound(repulsion_idx));

        // calculate slope points & minor correction to repulsion index
        vector<double> deriva;
        for (size_t i = repulsion_idx;
             i < repulsion_idx + 5 && i < energy.size() - 1; ++i) {
            double d = (energy[i + 1] - energy[i]) / __step_in_file;
            deriva.push_back(d);
        }

        // find up to 3 most negative derivatives and store index to slope
        set<size_t> slope;
        for (size_t i = 0; i < 3 && i < deriva.size(); ++i) {
            auto it = min_element(deriva.begin(), deriva.end(),
                                  [](double i, double j) { return i < j; });
            if (*it < 0) {  // derivative < 0
                int i0 = repulsion_idx + (it - deriva.begin());
                deriva.erase(it);
                slope.insert(i0);
                slope.insert(i0 + 1);
            }
        }

        try {
            if (slope.size() <= 1)
                throw InterpolationError("warning : slope not found in data");

            repulsion_idx = *slope.begin();  // correct repulsion_idx

            dbgmsg("atom1 = "
                   << idatm_type1 << " atom2 = " << idatm_type2 << " vdW_sum = "
                   << vdW_sum << " repulsion idx = " << repulsion_idx
                   << " step_non_bond = " << __step_non_bond
                   << " repulsion distance (below is forbidden area) = "
                   << __get_lower_bound(repulsion_idx) << " minimum distance = "
                   << __get_lower_bound(global_min_idx) << " begin slope idx = "
                   << *slope.begin() << " end slope idx = " << *slope.rbegin());

            // dbgmsg("energies before interpolations : " << endl
            //                                           << energy);

            const string idatm_type1 = help::idatm_unmask[atom_pair.first];
            const string idatm_type2 = help::idatm_unmask[atom_pair.second];

            vector<double> dataX, dataY;
            for (size_t i = repulsion_idx; i < energy.size();
                 ++i) {  // unset everything below this value
                dataX.push_back(__get_lower_bound(i));
                dataY.push_back(energy[i]);
            }

            Interpolation::BSplineFit BSFited(dataX, dataY);
            vector<double> potential = BSFited.interpolate_bspline(
                dataX.front(), dataX.back(), __step_non_bond);

            // Extrapolate extra points to fill the gap between the upperbound
            // and the cutoff
            // IE. If the cutoff is 6, the upperbound is 5.9, so 5.91 5.92, etc
            // need to be added
            // This is done here to avoid tainting the spline fitting
            const double potential_deriv_1 =
                (*(potential.rbegin() + 0) - *(potential.rbegin() + 1)) *
                __step_non_bond;
            const double potential_deriv_2 =
                (*(potential.rbegin() + 1) - *(potential.rbegin() + 2)) *
                __step_non_bond;
            const double potential_2_deriv =
                (potential_deriv_1 - potential_deriv_2) * __step_non_bond;

            for (size_t i = 0; i < (__step_in_file / __step_non_bond); ++i) {
                potential.push_back(potential.back() + (potential_deriv_1 +
                                                        i * potential_2_deriv));
            }

            // add repulsion term by fitting 1/x**12 function to slope points
            const double x1 = __get_lower_bound(*slope.begin());
            const double x2 = __get_lower_bound(*slope.rbegin());
            std::vector<double> r;
            std::vector<double> pot;
            int i = 0;
            for (double xi = dataX.front(); xi <= dataX.back();
                 xi += __step_non_bond) {
                if (xi >= x1 && xi <= x2) {
                    r.push_back(xi);
                    pot.push_back(potential[i]);
                }
                ++i;
            }

            dbgmsg("x1 = " << x1);
            dbgmsg("x2 = " << x2);

            // fit function to slope
            double coeffA, coeffB, WSSR;
            std::tie(coeffA, coeffB, WSSR) =
                fit_range_power_function_fast(r, pot);

            dbgmsg("atom1 = " << idatm_type1 << " atom2 = " << idatm_type2
                              << " coeffA = " << coeffA
                              << " coeffB = " << coeffB << " WSSR = " << WSSR);

            vector<double> repulsion;
            // Add starting point for x=0 (so it's its not NaN).
            repulsion.push_back(10 * coeffA / pow(__step_non_bond, 12) +
                                coeffB);

            // No longer loop over a double as this caused an extra point to be
            // calculated
            for (size_t i = 1; i < std::floor(__dist_cutoff / __step_non_bond) -
                                       potential.size() + 1;
                 ++i) {
                const double yi =
                    coeffA / pow(i * __step_non_bond, 12) + coeffB;
                repulsion.push_back(
                    std::isinf(yi)
                        ? 10 * coeffA / pow(i * (__step_non_bond + 1), 12) +
                              coeffB
                        : yi);
            }
#ifndef NDEBUG
            for (size_t i = 0; i < repulsion.size(); ++i)
                dbgmsg("i = " << i << " repulsion = " << repulsion[i]);
#endif
            // if the repulsion term comes under the potential do a linear
            // interpolation to get smooth joint
            int w = 0;
            while (!repulsion.empty() && repulsion.back() < potential.front()) {
                repulsion.pop_back();
                ++w;
            }

            dbgmsg("repulsion.size() = " << repulsion.size());

            const double rep_good = repulsion.back();
            const double k0 = (potential.front() - rep_good) / w;
            for (int k = 1; k <= w; ++k) {
                repulsion.push_back(rep_good + k0 * k);
            }

            // add repulsion term before bsplined potential
            potential.insert(potential.begin(), repulsion.begin(),
                             repulsion.end());

            __energies[atom_pair].assign(potential.begin(), potential.end());

            __derivatives[atom_pair] = Interpolation::derivative(
                __energies[atom_pair], __step_non_bond);

            // scale derivatives ONLY not energies and don't do it in kbforce
            // cause here is more efficient
            for (auto& dEdt : __derivatives[atom_pair]) {
                dEdt *= scale_non_bond;
            }

#ifndef NDEBUG
            for (size_t i = 0; i < potential.size(); ++i) {
                dbgmsg("interpolated "
                       << help::idatm_unmask[atom_pair.first] << " "
                       << help::idatm_unmask[atom_pair.second] << " "
                       << i * __step_non_bond << " pot = " << potential[i]);
            }
#endif
        }  // END of try
        catch (InterpolationError& e) {
            dbgmsg(e.what());
            const int n = (int)std::floor(__dist_cutoff / __step_non_bond) + 1;
            dbgmsg("n = " << n);
            __energies[atom_pair].assign(n, 0.0);
            __derivatives[atom_pair].assign(n, 0.0);
        }
    }
    dbgmsg("out of loop");
    log_benchmark << "Time to compile the objective function: "
                  << bench.seconds_from_start() << "\n";
    return *this;
}
}
}
