#include "score.hpp"
#include "interpolation.hpp"
#include "molib/molecule.hpp"
#include "helper/inout.hpp"
#include "helper/path.hpp"
#include "powerfit.hpp"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <functional>
#include <math.h>
#include <string>

namespace Score {
	ostream& operator<< (ostream& stream, const vector<double> &energy) {
		for (size_t i = 0; i < energy.size(); ++i) {
			stream << i << "=" << energy[i] << " ";
		}
		return stream;
	}

	ostream& operator<< (ostream& stream, const Score::AtomPairValues &energies) {
		for (auto &kv : energies) {
			auto &atom_pair = kv.first;
			auto &energy = kv.second;
			for (size_t i = 0; i < energy.size(); ++i) {
				stream << help::idatm_unmask[atom_pair.first] << " " << help::idatm_unmask[atom_pair.second] 
					<< " " << i << " " << energy[i] << endl;
			}
		}
		return stream;
	}

	ostream& operator<< (ostream& stream, const Score &score) {
		for (auto &kv : score.__energies_scoring) {
			auto &atom_pair = kv.first;
			auto &sene = score.__energies_scoring.at(atom_pair);
			auto &ene = score.__energies.at(atom_pair);
			auto &der = score.__derivatives.at(atom_pair);
			stream << "#objective function" << endl;
			for (size_t i = 0; i < ene.size(); ++i) {
				stream << help::idatm_unmask[atom_pair.first] << "\t"
					<< help::idatm_unmask[atom_pair.second] << "\t" 
					<< fixed << setprecision(3) << i * score.__step_non_bond << "\t" 
					<< fixed << setprecision(3) << ene[i] << "\t" 
					<< fixed << setprecision(3) << der[i] << endl;
			}
			stream << "#scoring function" << endl;
			for (size_t i = 0; i < sene.size(); ++i) {
				stream << help::idatm_unmask[atom_pair.first] << "\t"
					<< help::idatm_unmask[atom_pair.second] << "\t" 
					<< fixed << setprecision(3) << score.__get_lower_bound(i) << "\t" 
					<< fixed << setprecision(3) << sene[i] << endl;
			}
		}
		return stream;
	}

	Score& Score::output_objective_function(const string &obj_dir) {
		for (auto &kv : __energies) {
			auto &atom_pair = kv.first;
			auto &ene = __energies.at(atom_pair);
			//~ auto &der = __derivatives.at(atom_pair);

			stringstream ss;

			const string idatm_type1 = help::idatm_unmask[atom_pair.first];
			const string idatm_type2 = help::idatm_unmask[atom_pair.second];

			for (size_t i = 0; i < ene.size(); ++i) {
				ss << setprecision(8) << ene[i] << endl;
			}
			const string &filename = idatm_type1 + "_" + idatm_type2 + ".txt";
			
			Inout::file_open_put_stream(Path::join(Path::join(obj_dir, std::to_string(__step_non_bond)), filename), ss);
		}
		return *this;
	}

        Score& Score::parse_objective_function(const string &obj_dir, const double scale_non_bond, const size_t max_step) {

                log_step << "Parsing obj function from disk" << endl;

                boost::filesystem::path path_to_objective_function(obj_dir);
                std::string subdir = std::to_string(__step_non_bond);

                // Some systems add trailing zeros to doubles, let's remove them
                while ( ! boost::filesystem::exists(path_to_objective_function / subdir) && subdir.back() == '0' ) {
                        subdir.erase(subdir.end() - 1);
                }

                path_to_objective_function /= subdir;

                if ( ! boost::filesystem::exists(path_to_objective_function) ) {
                        throw Error( "Objective function not found! check 'obj_dir' and 'step'");
                }

                for (auto &atom_pair : __avail_prot_lig) {

                        const string idatm_type1 = help::idatm_unmask[atom_pair.first];
                        const string idatm_type2 = help::idatm_unmask[atom_pair.second];

                        dbgmsg("parsing objective function for " << idatm_type1 << " and " << idatm_type2);

                        const string &filename = idatm_type1 + "_" + idatm_type2 + ".txt";

                        // Check if the file exists on disk. Since we know that the directory exists on disk,
                        // we can skip this pair as it is not part of the scoring function (ie not in the CSD).
                        if ( ! Inout::file_size( (path_to_objective_function / filename).string() )) {
                                continue;
                        }

                        vector<string> contents;
                        Inout::read_file( (path_to_objective_function / filename).string(), contents);

                        for (auto &line : contents) {
                                stringstream ss(line);
                                string str1;
                                ss >> str1;
                                __energies[atom_pair].push_back(stod(str1));
                        }

                        if ( __energies[atom_pair].size() < max_step ) {
                                int energies_diff = max_step - __energies[atom_pair].size();

                                if (energies_diff < 0) {
                                        log_warning << "Large atom pair: " << __energies[atom_pair].size()
                                                    << " " << idatm_type1 << " " << idatm_type2 << endl;
                                        __energies[atom_pair].resize(max_step);
                                }

                                for ( int i = 0; i < energies_diff; ++i ) {
                                        __energies[atom_pair].push_back(0.0);
                                }
                        }

                        __derivatives[atom_pair] = Interpolation::derivative(__energies[atom_pair], __step_non_bond);

                        // scale derivatives ONLY not energies and don't do it in kbforce cause here is more efficient
                        for (auto &dEdt : __derivatives[atom_pair]) {
                                dEdt *= scale_non_bond;
                        }

                }
                dbgmsg("parsed objective function");
                return *this;
        }

        Array1d<double> Score::compute_energy(const Molib::Atom::Grid &gridrec, const Geom3D::Coordinate &crd, const set<int> &ligand_atom_types) const {
                dbgmsg("computing energy");
                Array1d<double> energy_sum(*ligand_atom_types.rbegin() + 1);
                for (auto &patom : gridrec.get_neighbors(crd, __dist_cutoff)) {
                        const double dist = patom->crd().distance(crd);
                        const auto &atom_1 = patom->idatm_type();
                        const int index = __get_index(dist);
                        for (auto &l : ligand_atom_types) {
                                auto atom_pair = minmax(atom_1, l);
#ifndef NDEBUG
                                if (!__energies_scoring.count(atom_pair))
                                        throw Error("undefined atom_pair in __energies_scoring");
                                if (static_cast<size_t> (index) >= __energies_scoring.at(atom_pair).size())
                                        throw Error("undefined index in __energies_scoring");
                                dbgmsg("atom pairs = " << help::idatm_unmask[atom_pair.first] << " " 
                                        << help::idatm_unmask[atom_pair.second] << " index = " << index
                                        << " dist = " << dist << " __energies_scoring.at(atom_pair).size() = " 
                                        << __energies_scoring.at(atom_pair).size() << " partial energy = "
                                        << __energies_scoring.at(atom_pair).at(index));
#endif
                                if ( static_cast<size_t>(index) >= __energies_scoring.at(atom_pair).size()) {
                                        log_warning << "An index from get_neighbors is greater than cutoff by step_in_file ("
                                                    << index << " and " << __dist_cutoff / __step_in_file << " respectively).\n"
                                                    << " -- Original distance is " << dist << "\n"
                                                    << " -- Setting point's score to zero." << endl;
                                        continue;
                                }
                                energy_sum.data[l] += __energies_scoring.at(atom_pair).at(index);
                        }
                }
                dbgmsg("out of compute energy energy_sum = " << energy_sum);
                return energy_sum;
        }

        Score& Score::define_composition(const set<int> &receptor_idatm_types, const set<int> &ligand_idatm_types) {
                set<int> idatm_types;
                for (auto &key : receptor_idatm_types) idatm_types.insert(key);
                for (auto &key : ligand_idatm_types) idatm_types.insert(key);

                dbgmsg("idatm_types.size() = " << idatm_types.size());
                for (auto &prot_key : idatm_types) {
                        for (auto &lig_key : idatm_types) {
                                __prot_lig_pairs.insert(minmax(prot_key, lig_key));
                        }
                }

                __avail_prot_lig = __prot_lig_pairs;

                if (__comp == "complete") {
                        int sz = help::idatm_mask.size();
                        for (int i = 0; i < sz; ++i) {
                                for (int j = i; j < sz; ++j) {
                                        __prot_lig_pairs.insert({i, j});
                                }
                        }
                }

                dbgmsg("__prot_lig_pairs.size() = " << __prot_lig_pairs.size());
#ifndef NDEBUG
                for (auto &i : __prot_lig_pairs) {
                        dbgmsg("pairs: " << help::idatm_unmask[i.first] << " " 
                                << help::idatm_unmask[i.second]); 
                }
#endif
                return *this;
        }

	Score& Score::process_distributions_file(const string &distributions_file) {
		Benchmark bench;
		log_step << "processing combined histogram ...\n";
		vector<string> distributions_file_raw;
		Inout::read_file(distributions_file, distributions_file_raw);
		const bool rad_or_raw(__rad_or_raw == "normalized_frequency");

		for (string &line : distributions_file_raw) {
			stringstream ss(line); // dist_file is simply too big to use boost
			string atom_1,atom_2;
			double lower_bound, upper_bound, quantity;
			ss>>atom_1>>atom_2>>lower_bound>>upper_bound>>quantity;
			
			// set step size
			if (__step_in_file < 0) {
				__step_in_file = upper_bound - lower_bound;
				dbgmsg("step_in_file = " << __step_in_file);
				__bin_range_sum.resize(__get_index(__dist_cutoff) + 1, 0);
			}
				
			
			if (upper_bound <= __dist_cutoff) {
				pair_of_ints atom_pair = minmax(help::idatm_mask.at(atom_1), help::idatm_mask.at(atom_2));
				if (__prot_lig_pairs.count(atom_pair)) {
					double shell_volume = (rad_or_raw ? 1.0 : 4*M_PI*pow(upper_bound, 3)/3 - 4*M_PI*pow(lower_bound, 3)/3);
					__gij_of_r_numerator[atom_pair].push_back(quantity / shell_volume);
					__sum_gij_of_r_numerator[atom_pair] += quantity / shell_volume;
					// JANEZ : next two are for cumulative scoring function (compile_cumulative_scoring_function)
					__bin_range_sum[__get_index(lower_bound)] += quantity / shell_volume;
					__total_quantity += quantity / shell_volume;
					//~ dbgmsg(" " <<  atom_1 << " " <<  atom_2 << " " << lower_bound 
						//~ << " " <<  upper_bound << " " <<  quantity);
				}
			}	
		}
		// JANEZ : next part only needed for compile_mean_scoring_function
		__gij_of_r_bin_range_sum.resize(__get_index(__dist_cutoff) + 1, 0);
		for (auto &el1 : __gij_of_r_numerator) {
			const pair_of_ints &atom_pair = el1.first;
			if (__sum_gij_of_r_numerator[atom_pair] > 0) {
				const vector<double> &gij_of_r_vals = el1.second;
				for (size_t i = 0; i < gij_of_r_vals.size(); ++i) {
					const double &gij_of_r_value = gij_of_r_vals[i];
					__gij_of_r_bin_range_sum[i] += gij_of_r_value / __sum_gij_of_r_numerator[atom_pair];
					dbgmsg("__gij_of_r_bin_range_sum[" <<  __get_lower_bound(i) 
						<< "]= " <<  __gij_of_r_bin_range_sum[i]);
				}
			}
		}
		log_benchmark << "time to process distributions file " << bench.seconds_from_start() 
			<< " wallclock seconds" << "\n";
		return *this;
	}

	Score& Score::compile_objective_function(const double scale_non_bond) {
		log_step << "Compiling objective function for minimization...\n";
                Benchmark bench;
		auto energy_function = __ref_state == "mean" ? 
			mem_fn(&Score::__energy_mean) : mem_fn(&Score::__energy_cumulative);
		for (auto &el1 : __gij_of_r_numerator) {
			const pair_of_ints &atom_pair = el1.first;

                        // Don't bother making an objective function
                        // unless we're going to use it
                        if (!__avail_prot_lig.count(atom_pair))
                                continue;
                        
			const vector<double> &gij_of_r_vals = el1.second;
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
				const double &gij_of_r_numerator = gij_of_r_vals[i];
				dbgmsg("lower bound for atom_pair " << idatm_type1
					<< " " << idatm_type2 
					<< " " << lower_bound << " gij_of_r_numerator = "
					<< gij_of_r_numerator << " __sum_gij_of_r_numerator[atom_pair] = "
					<< __sum_gij_of_r_numerator[atom_pair]);

				if (__sum_gij_of_r_numerator[atom_pair] < __eps) {
					energy[i] = 0;
				}
				else if (gij_of_r_numerator < __eps && (i + 1 < start_idx)) {
					energy[i] = 5.0;
				}
				else if (gij_of_r_numerator < __eps && (i + 1 >= start_idx)) {
					energy[i] = 0;
				}
				else {
					energy[i] = energy_function(*this, atom_pair, lower_bound);
				}
			}

			dbgmsg("raw energies before interpolations : " << endl << energy);

			// correct for outliers only in the TRUE potential region of the potential
			for (size_t i = start_idx; i < energy.size() - 2; ++i) {
				if (energy[i] != 0 && energy[i] != 5) {
					size_t j = i + 1;
					while (j < i + 3 && j < energy.size() && (energy[j] == 0 || energy[j] == 5)) { ++j; }
					if (j > i + 1 && j < energy.size()) {
						const double k0 = (energy[j] - energy[i]) / (j - i);
						for (size_t k = 1; k < j - i; ++k) {
							energy[i + k] = energy[i] + k0 * k;
						}
					}
				}
			}

			// locate global minimum and repulsion index in interval [vdW_sum - 0.6, vdW_sum + 1.0]
			double global_min = HUGE_VAL;
			size_t global_min_idx = start_idx;
			size_t repulsion_idx = global_min_idx;

			try { 
				repulsion_idx = __get_index(help::repulsion_idx.at(make_pair(idatm_type1, idatm_type2)));
			} catch (out_of_range&) {
				try {
					repulsion_idx = __get_index(help::repulsion_idx.at(make_pair(idatm_type2, idatm_type1))); 
				} catch (out_of_range&) {

					dbgmsg("de-novo calculation of repulsion idx");
					
					for (size_t i = start_idx; i < end_idx && i < energy.size(); ++i) {
						if (energy[i] != -HUGE_VAL) {
							if (energy[i] < global_min) {
								global_min = energy[i];
								global_min_idx = i;
								
								// minimum has to have steep downward slope on the left side
								// leave the slope intact and unset everything left of the slope start point
								repulsion_idx = i;
								while (repulsion_idx > 0 && energy[repulsion_idx] != -HUGE_VAL 
									&& ((energy[repulsion_idx] - energy[repulsion_idx - 1]) / __step_in_file) < 0.75) --repulsion_idx;
		
							}
						}
					}

				}
			}

			dbgmsg("repulsion idx (before correction) = " << repulsion_idx
				<< " repulsion distance (below is forbidden area) = " << __get_lower_bound(repulsion_idx));
			
			// calculate slope points & minor correction to repulsion index
			vector<double> deriva;
			for (size_t i = repulsion_idx; i < repulsion_idx + 5 && i < energy.size() - 1; ++i) {
				double d = (energy[i + 1] - energy[i]) / __step_in_file;
				deriva.push_back(d);
			}
			
			// find up to 3 most negative derivatives and store index to slope
			set<size_t> slope;
			for (size_t i = 0; i < 3 && i < deriva.size(); ++i) {
				auto it = min_element(deriva.begin(), deriva.end(), [](double i, double j) { return i < j; });
				if (*it < 0) { // derivative < 0
					int i0 = repulsion_idx + (it - deriva.begin());
					deriva.erase(it);
					slope.insert(i0);
					slope.insert(i0 + 1);
				}
			}
			
			try {

				if (slope.size() <= 1) throw InterpolationError("warning : slope not found in data");

				repulsion_idx = *slope.begin(); // correct repulsion_idx
			
				dbgmsg("atom1 = " << idatm_type1
					<< " atom2 = " << idatm_type2
					<< " vdW_sum = " << vdW_sum
					<< " repulsion idx = " << repulsion_idx
					<< " step_non_bond = " << __step_non_bond
					<< " repulsion distance (below is forbidden area) = " << __get_lower_bound(repulsion_idx)
					<< " minimum distance = " << __get_lower_bound(global_min_idx)
					<< " begin slope idx = " << *slope.begin()
					<< " end slope idx = " << *slope.rbegin()
				);
	
                                dbgmsg("energies before interpolations : " << endl << energy);
    
                                const string idatm_type1 = help::idatm_unmask[atom_pair.first];
                                const string idatm_type2 = help::idatm_unmask[atom_pair.second];

                                vector<double> dataX, dataY;
                                for (size_t i = repulsion_idx; i < energy.size(); ++i) { // unset everything below this value
                                        dataX.push_back(__get_lower_bound(i));
                                        dataY.push_back(energy[i]);
                                }

                                Interpolation::BSplineFit BSFited(dataX, dataY);
                                vector<double> potential = BSFited.interpolate_bspline(dataX.front(), dataX.back(), __step_non_bond);

                                // Extrapolate extra points to fill the gap between the upperbound and the cutoff
                                // IE. If the cutoff is 6, the upperbound is 5.9, so 5.91 5.92, etc need to be added
                                // This is done here to avoid tainting the spline fitting
                                const double potential_deriv_1 = (*(potential.rbegin() + 0) - *(potential.rbegin() + 1)) * __step_non_bond;
                                const double potential_deriv_2 = (*(potential.rbegin() + 1) - *(potential.rbegin() + 2)) * __step_non_bond;
                                const double potential_2_deriv = (potential_deriv_1 - potential_deriv_2) * __step_non_bond;

                                for ( size_t i = 0; i < (__step_in_file / __step_non_bond); ++i) {
                                        potential.push_back(potential.back() + (potential_deriv_1 + i*potential_2_deriv ));
                                }

                                // add repulsion term by fitting 1/x**12 function to slope points
                                const double x1 = __get_lower_bound(*slope.begin());
                                const double x2 = __get_lower_bound(*slope.rbegin());
                                std::vector<double> r;
                                std::vector<double> pot;
                                int i = 0;
                                for (double xi = dataX.front(); xi <= dataX.back(); xi += __step_non_bond) {
                                        if (xi >= x1 && xi <=x2) {
                                                r.push_back(xi);
                                                pot.push_back(potential[i]);
                                        }
                                        ++i;
                                }

                                dbgmsg("x1 = " << x1);
                                dbgmsg("x2 = " << x2);

                                // fit function to slope
                                double coeffA, coeffB, WSSR;
                                std::tie(coeffA, coeffB, WSSR) = fit_range_power_function_fast(r, pot);

                                dbgmsg("atom1 = " << idatm_type1
                                        << " atom2 = " << idatm_type2
                                        << " coeffA = " << coeffA
                                        << " coeffB = " << coeffB
                                        << " WSSR = " << WSSR);


                                vector<double> repulsion;
                                // Add starting point for x=0 (so it's its not NaN).
                                repulsion.push_back(10 * coeffA / pow(__step_non_bond , 12) + coeffB);

                                // No longer loop over a double as this caused an extra point to be calculated
                                for (size_t i = 1; i < std::floor( __dist_cutoff / __step_non_bond) - potential.size() + 1; ++i) {
                                        const double yi = coeffA / pow(i * __step_non_bond, 12) + coeffB;
                                        repulsion.push_back(std::isinf(yi) ? 10 * coeffA / pow(i * (__step_non_bond+1) , 12) + coeffB : yi);
                                }
#ifndef NDEBUG
				for (size_t i = 0; i < repulsion.size(); ++i)
					dbgmsg("i = " << i << " repulsion = " << repulsion[i]);
#endif
				// if the repulsion term comes under the potential do a linear interpolation to get smooth joint
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
				potential.insert(potential.begin(), repulsion.begin(), repulsion.end());

                                __energies[atom_pair].assign(potential.begin(), potential.end());
                                
                                __derivatives[atom_pair] = Interpolation::derivative(__energies[atom_pair], __step_non_bond);

                                // scale derivatives ONLY not energies and don't do it in kbforce cause here is more efficient
                                for (auto &dEdt : __derivatives[atom_pair]) {
                                        dEdt *= scale_non_bond;
                                }

#ifndef NDEBUG
				for (size_t i = 0; i < potential.size(); ++i) {
					dbgmsg("interpolated " << help::idatm_unmask[atom_pair.first] 
						<< " " << help::idatm_unmask[atom_pair.second]
						<< " " << i * __step_non_bond
						<< " pot = " << potential[i]);
				}
#endif
			} // END of try
			catch (InterpolationError &e) {
				dbgmsg(e.what());
				const int n = (int) std::floor( __dist_cutoff / __step_non_bond) + 1;
				dbgmsg("n = " << n);
				__energies[atom_pair].assign(n, 0.0);
                                __derivatives[atom_pair].assign(n,0.0);
			}
		}
		dbgmsg("out of loop");
                log_benchmark << "Time to compile the objective function: " << bench.seconds_from_start() << "\n";
		return *this;
	}
	
	Score& Score::compile_scoring_function() {
		log_step << "Compiling scoring function...\n";
		auto energy_function = __ref_state == "mean" ? 
			mem_fn(&Score::__energy_mean) : mem_fn(&Score::__energy_cumulative);
		for (auto &el1 : __gij_of_r_numerator) {
			const pair_of_ints &atom_pair = el1.first;
			const vector<double> &gij_of_r_vals = el1.second;
			const double w1 = help::vdw_radius[atom_pair.first];
			const double w2 = help::vdw_radius[atom_pair.second];
			const double vdW_sum = ((w1 > 0 && w2 > 0) ? w1 + w2 : 4.500);
			const size_t repulsion_idx = __get_index(vdW_sum - 0.6);
			const string idatm_type1 = help::idatm_unmask[atom_pair.first];
			const string idatm_type2 = help::idatm_unmask[atom_pair.second];
			dbgmsg("atom1 = " << idatm_type1
				<< " atom2 = " << idatm_type2
				<< " vdW_sum = " << vdW_sum
				<< " repulsion index (below is forbidden area) = " << repulsion_idx);
			vector<double> energy(gij_of_r_vals.size(), -HUGE_VAL);
			for (size_t i = 0; i < gij_of_r_vals.size(); ++i) {
				const double lower_bound = __get_lower_bound(i);
				const double &gij_of_r_numerator = gij_of_r_vals[i];
				dbgmsg("lower bound for atom_pair " << idatm_type1
					<< " " << idatm_type2 
					<< " " << lower_bound << " gij_of_r_numerator = "
					<< gij_of_r_numerator << " __sum_gij_of_r_numerator[atom_pair] = "
					<< __sum_gij_of_r_numerator[atom_pair]);

				if (__sum_gij_of_r_numerator[atom_pair] < __eps) {
					energy[i] = 0;
				}
				else if (gij_of_r_numerator < __eps && (i + 1 < repulsion_idx)) {
					energy[i] = 5.0;
				}
				else if (gij_of_r_numerator < __eps && (i + 1 >= repulsion_idx)) {
					energy[i] = 0;
				}
				else {
					energy[i] = energy_function(*this, atom_pair, lower_bound);
				}
			}

			dbgmsg("raw energies for scoring : " << endl << energy);


			if (idatm_type1 == "H" || idatm_type1 == "HC"
				|| idatm_type2 == "H" || idatm_type2 == "HC") {
				energy.assign(energy.size(), 0);
			}
			__energies_scoring[atom_pair] = energy;
		}
		dbgmsg("out of loop");
		return *this;
	}
	
	double Score::__energy_mean(const pair_of_ints &atom_pair, const double &lower_bound) {
		const int idx = __get_index(lower_bound);
#ifndef NDEBUG
		if (!__gij_of_r_numerator.count(atom_pair))
			throw Error("die : undefined __gij_of_r_numerator");
		if (static_cast<size_t>(idx) >= __gij_of_r_numerator.at(atom_pair).size())
			throw Error("die : undefined __gij_of_r_numerator");
		if (static_cast<size_t>(idx) >= __gij_of_r_bin_range_sum.size())
			throw Error("die : undefined __gij_of_r_bin_range_sum");
#endif
		double gij_of_r = __gij_of_r_numerator[atom_pair][idx] / __sum_gij_of_r_numerator[atom_pair];
		dbgmsg("gij_of_r = " << gij_of_r);
		double denominator = __gij_of_r_bin_range_sum[idx] / __prot_lig_pairs.size();
		dbgmsg("denominator = " << denominator);
		double ratio = gij_of_r / denominator;
		return -log(ratio);
	}
	
	double Score::__energy_cumulative(const pair_of_ints &atom_pair, const double &lower_bound) {
		const int idx = __get_index(lower_bound);
#ifndef NDEBUG
		if (!__gij_of_r_numerator.count(atom_pair))
			throw Error("die : undefined __gij_of_r_numerator");
		if (static_cast<size_t>(idx) >= __gij_of_r_numerator.at(atom_pair).size())
			throw Error("die : undefined __gij_of_r_numerator");
		if (static_cast<size_t>(idx) >= __bin_range_sum.size())
			throw Error("die : undefined __bin_range_sum");
#endif
		double numerator = __gij_of_r_numerator[atom_pair][idx] / __sum_gij_of_r_numerator[atom_pair];
		double denominator = __bin_range_sum[idx] / __total_quantity;
		double ratio = numerator / denominator;
		return -log(ratio);
	}
	
	double Score::non_bonded_energy(const Molib::Atom::Grid &gridrec, const Molib::Molecule &ligand) const {
		return this->non_bonded_energy(gridrec, ligand.get_atoms(), ligand.get_crds());
	}

	double Score::non_bonded_energy(const Molib::Atom::Grid &gridrec, const Molib::Atom::Vec &atoms, const Geom3D::Point::Vec &crds) const {
		double energy_sum = 0.0;
		for (size_t i = 0; i < atoms.size(); ++i) {
			const Molib::Atom &atom2 = *atoms[i];
			const Geom3D::Coordinate &atom2_crd = crds[i];
			const auto &atom_2 = atom2.idatm_type();
			dbgmsg("ligand atom = " << atom2.atom_number() << " crd= " << atom2_crd.pdb());
			for (auto &atom1 : gridrec.get_neighbors(atom2_crd, __dist_cutoff)) {
				const double dist = atom1->crd().distance(atom2_crd);
				dbgmsg("dist = " << setprecision(12) << dist);
				dbgmsg("dist_sq = " << setprecision(12) << atom1->crd().distance_sq(atom2_crd));
				const auto &atom_1 = atom1->idatm_type();
				const pair_of_ints atom_pair = minmax(atom_1, atom_2);
				const int idx = __get_index(dist);
				energy_sum += __energies_scoring.at(atom_pair).at(idx);
#ifndef NDEBUG				
				dbgmsg("ligand atom = " << atom2.atom_number() 
					<< " crd= " << atom2_crd.pdb() << "protein atom = " << atom1->atom_number() 
					<< " crd= " << atom1->crd().pdb() << " dist= " << dist 
					<< " atom_pair={" << help::idatm_unmask[atom_pair.first] << "," << help::idatm_unmask[atom_pair.second] << "}" 
					<< " lower_bound= " << __get_lower_bound(idx) << " energy_sum=" << energy_sum);
#endif
			}
		}
		dbgmsg("exiting non_bonded_energy");
		return energy_sum;
	}

};
