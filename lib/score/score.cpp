#include "score.hpp"
#include "molib/molecule.hpp"
#include "helper/inout.hpp"
#include "helper/path.hpp"
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
				const size_t idx = __get_index(dist);

				// The index function can produce values that are too large due to the +0.0000000001
				// We ignore these values because they are technically > the cutoff, if only by a little
				if (idx < __energies_scoring.at(atom_pair).size()) {
					energy_sum += __energies_scoring.at(atom_pair).at(idx);
				} else { 
					log_warning << "Close call: "
						    << setprecision(12) << dist << " "
						    << atom1->crd() << "\t" << atom2_crd << " "
						    << "difference: " << __energies_scoring.at(atom_pair).back()
						    << endl;
				}
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
