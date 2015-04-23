#include "score.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/inout.hpp"
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace Molib {
	ostream& operator<< (ostream& stream, const Score::M0 &energy) {
		for (int i = 0; i < energy.size(); ++i) {
			stream << "i = " << i << " value = " << energy[i] << endl;
		}
		return stream;
	}

	double Score::__distances_and_scores(const Geom3D::Coordinate &crd2, int ligand_atom_type) {
		double energy_sum = 0.0;
		Atom at(crd2);
		for (auto &patom : __gridrec.get_neighbors(at, __dist_cutoff)) {
			const double dist = patom->crd().distance(crd2);
			const auto &atom_1 = patom->idatm_type();
			auto atom_pair = minmax(atom_1, ligand_atom_type);
			const int index = __get_index(dist);
#ifndef NDEBUG
			dbgmsg("atom pairs = " << atom_pair.first << " " 
				<< atom_pair.second << " index = " << index
				<< " dist = " << dist << " __energies.at(atom_pair).size() = " 
				<< __energies.at(atom_pair).size());
			if (!__energies.count(atom_pair))
				throw Error("undefined atom_pair in __energies");
			if (index >= __energies.at(atom_pair).size())
				throw Error("undefined index in __energies");
#endif
			energy_sum += __energies[atom_pair][index];
		}
		dbgmsg("out of distances and scores");
		return energy_sum;
	}
	void Score::__define_composition(const set<int> &receptor_idatm_types, const set<int> &ligand_idatm_types) {
		// JANEZ: num_pairs is now __prot_lig_pairs.size() !!!
		set<int> idatm_types;
		for (auto &key : receptor_idatm_types) idatm_types.insert(key);
		for (auto &key : ligand_idatm_types) idatm_types.insert(key);
		if (__comp == "reduced") {
			for (auto &prot_key : idatm_types) {
				for (auto &lig_key : idatm_types) {
					__prot_lig_pairs.insert(minmax(prot_key, lig_key));
				}
			}
		} else if (__comp == "complete") {
			int sz = sizeof(help::idatm_unmask) / sizeof(help::idatm_unmask[0]);
			for (int i = 0; i < sz; ++i) {
				for (int j = i + 1; j < sz; ++j) {
					__prot_lig_pairs.insert({i, j});
				}
			}
		}
#ifndef NDEBUG
		for (auto &i : __prot_lig_pairs) {
			dbgmsg("pairs: " << help::idatm_unmask[i.first] << " " 
				<< help::idatm_unmask[i.second]); 
		}
#endif
	}
	void Score::__process_distributions_file() {
		Benchmark::reset();
		cout << "processing combined histogram ...\n";
		vector<string> distributions_file_raw;
		inout::Inout::read_file(__distributions_file, distributions_file_raw);
		const bool rad_or_raw(__rad_or_raw == "normalized_frequency");
		__bin_range_sum.resize(__get_index(__dist_cutoff) + 1, 0);
		for (string &line : distributions_file_raw) {
			stringstream ss(line); // dist_file is simply too big to use boost
			string atom_1,atom_2;
			double lower_bound, upper_bound, quantity;
			ss>>atom_1>>atom_2>>lower_bound>>upper_bound>>quantity;
			if (upper_bound <= __dist_cutoff) {
				pair_of_ints atom_pair = minmax(help::idatm_mask.at(atom_1), help::idatm_mask.at(atom_2));
				if (__prot_lig_pairs.count(atom_pair)) {
					double shell_volume = (rad_or_raw ? 1.0 : 4*M_PI*pow(upper_bound, 3)/3 - 4*M_PI*pow(lower_bound, 3)/3);
					__gij_of_r_numerator[atom_pair].push_back(quantity / shell_volume);
					__sum_gij_of_r_numerator[atom_pair] += quantity / shell_volume;
					// JANEZ : next two are for cumulative scoring function (compile_cumulative_scoring_function)
					__bin_range_sum[__get_index(lower_bound)] += quantity / shell_volume;
					__total_quantity += quantity / shell_volume;
					dbgmsg(" " <<  atom_1 << " " <<  atom_2 << " " << lower_bound 
						<< " " <<  upper_bound << " " <<  quantity);
				}
			}	
		}
		// JANEZ : next part only needed for compile_mean_scoring_function
		__gij_of_r_bin_range_sum.resize(__get_index(__dist_cutoff) + 1, 0);
		for (auto &el1 : __gij_of_r_numerator) {
			const pair_of_ints &atom_pair = el1.first;
			if (__sum_gij_of_r_numerator[atom_pair] > 0) {
				const M0 &gij_of_r_vals = el1.second;
				for (int i = 0; i < gij_of_r_vals.size(); ++i) {
					const double &gij_of_r_value = gij_of_r_vals[i];
					__gij_of_r_bin_range_sum[i] += gij_of_r_value / __sum_gij_of_r_numerator[atom_pair];
					dbgmsg("__gij_of_r_bin_range_sum[" <<  __get_lower_bound(i) 
						<< "]= " <<  __gij_of_r_bin_range_sum[i]);
				}
			}
		}
		cout << "time to process distributions file " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
	}
	Score::Inter Score::__interpolate(const M0 &energy) {
		M0 x, y;
		Inter inter;
		for (int i = 0; i < energy.size(); ++i) {
			if (energy[i] != -HUGE_VAL) {
				x.push_back(__get_lower_bound(i));
				y.push_back(energy[i]);
			}
#ifndef NDEBUG
			else {
				dbgmsg("NOT FOR INTERPOLATION i = " << i 
					<< " energy = " << energy[i]);
			}
#endif
		}
		if (x.size() < 10) { // ?#! this is crappy, not enough points for interpolation
			dbgmsg("not enough points for interpolation, just zeroing everything!");
			for (int i = 0; i < energy.size(); ++i) {
				inter.potential.push_back(0.0);
			}
			for (int i = 0; i < inter.potential.size() - 1; ++i) {
				inter.derivative.push_back(0.0);
			}
		} else {
			gsl_interp_accel *acc = gsl_interp_accel_alloc();
			gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
			gsl_spline_init(spline, &x[0], &y[0], x.size());
			dbgmsg("min_bound = " << x.front() << " max_bound = " << x.back());
			for (int i = 0; i < energy.size(); ++i) {
				inter.potential.push_back(gsl_spline_eval(spline, __get_lower_bound(i), acc));
			}
			for (int i = 0; i < inter.potential.size() - 1; ++i) {
				// gsl_spline_eval_deriv gave WRONG derivatives ???
				inter.derivative.push_back((inter.potential[i + 1] - inter.potential[i]) / __step_non_bond);
			}
			inter.derivative.push_back(0); // last
			gsl_spline_free(spline);
			gsl_interp_accel_free(acc);
		}
		return inter;
	}

	void Score::__compile_scoring_function() {
		cout << "compiling scoring function ...\n";
		dbgmsg("her i am");
		auto energy_function = __ref_state == "mean" ? 
			mem_fn(&Score::__energy_mean) : mem_fn(&Score::__energy_cumulative);
		for (auto &el1 : __gij_of_r_numerator) {
			const pair_of_ints &atom_pair = el1.first;
			const M0 &gij_of_r_vals = el1.second;
			const double w1 = help::vdw_radius[atom_pair.first];
			const double w2 = help::vdw_radius[atom_pair.second];
			const double vdW_sum = ((w1 > 0 && w2 > 0) ? w1 + w2 : 4.500);
			const int repulsion_idx = __get_index(vdW_sum - 0.6);
			const string idatm_type1 = help::idatm_unmask[atom_pair.first];
			const string idatm_type2 = help::idatm_unmask[atom_pair.second];
			dbgmsg("atom1 = " << idatm_type1
				<< " atom2 = " << idatm_type2
				<< " vdW_sum = " << vdW_sum
				<< " repulsion index (below is forbidden area) = " << repulsion_idx);
			M0 energy(gij_of_r_vals.size(), -HUGE_VAL);
			for (int i = 0; i < gij_of_r_vals.size(); ++i) {
				const double lower_bound = __get_lower_bound(i);
				const double &gij_of_r_numerator = gij_of_r_vals[i];
				dbgmsg("lower bound for atom_pair " << idatm_type1
					<< " " << idatm_type2 
					<< " " << lower_bound << " gij_of_r_numerator = "
					<< gij_of_r_numerator << " __sum_gij_of_r_numerator[atom_pair] = "
					<< __sum_gij_of_r_numerator[atom_pair]);
				if (gij_of_r_numerator >= __eps)
					energy[i] = energy_function(*this, atom_pair, lower_bound);
			}
			//~ dbgmsg("maximum is at " << __get_lower_bound(max_index) 
				//~ << " and is " << max_energy);
			dbgmsg("raw energies before interpolations : " << endl << energy);
			// correct for outliers
			for (int i = 0; i < energy.size() - 2; ++i) {
				if (energy[i] != -HUGE_VAL && energy[i + 1] == -HUGE_VAL
					&& energy[i + 2] != -HUGE_VAL) {
					energy[i + 1] = (energy[i] + energy[i + 2]) / 2;
				}
				if (energy[i] == HUGE_VAL) energy[i] = -HUGE_VAL;
			}
			//~ if (max_index == -1 || __get_lower_bound(max_index) > 4.5 || __get_lower_bound(max_index) < 2.0) { // there are no data for this atom type
				//~ for (int i = 0; i < gij_of_r_vals.size(); ++i)
					//~ energy[i] = 0;
			//~ } else {
			//~ if (max_index == -1) { // there are no data for this atom type
				//~ for (int i = 0; i < gij_of_r_vals.size(); ++i)
					//~ energy[i] = 0;
			//~ } else {
				//~ for (int i = 0; i <= max_index; ++i) // unset everything below this value
					//~ energy[i] = -HUGE_VAL;
				//~ energy.front() = 280000; // set the lower bound energy to arbitrary high value
				//~ if (energy.back() == -HUGE_VAL) energy.back() = 0; // take care about the non-assigned values
			//~ }
			//~ if (max_index == -1) { // there are no data for this atom type
				//~ for (int i = 0; i < gij_of_r_vals.size(); ++i)
					//~ energy[i] = 0;
			//~ } else {
			for (int i = 0; i <= repulsion_idx; ++i) // unset everything below this value
				energy[i] = -HUGE_VAL;
			energy.front() = 280000; // set the lower bound energy to arbitrary high value
			if (energy.back() == -HUGE_VAL) energy.back() = 0; // take care about the non-assigned values
			if (idatm_type1 == "H" || idatm_type1 == "HC"
				|| idatm_type2 == "H" || idatm_type2 == "HC") {
				energy.assign(energy.size(), 0);
			}
			dbgmsg("energies before interpolations : " << endl << energy);
			auto inter = __interpolate(energy); // interpolate data points to get a smooth function
			__energies[atom_pair].assign(inter.potential.begin(), inter.potential.end());
			__derivatives[atom_pair].assign(inter.derivative.begin(), inter.derivative.end());
#ifndef NDEBUG
			for (int i = 0; i < inter.potential.size(); ++i) {
				dbgmsg("inter.potential " << help::idatm_unmask[atom_pair.first] 
					<< " " << help::idatm_unmask[atom_pair.second]
					<< " " << __get_lower_bound(i) << " " << inter.potential[i]); 
			}
			for (int i = 0; i < inter.derivative.size(); ++i) {
				dbgmsg("inter.derivative " << help::idatm_unmask[atom_pair.first] 
					<< " " << help::idatm_unmask[atom_pair.second]
					<< " " << __get_lower_bound(i) << " " << inter.derivative[i]); 
			}
#endif
		}
		dbgmsg("out of loop");
	}
	double Score::__energy_mean(const pair_of_ints &atom_pair, const double &lower_bound) {
		const int idx = __get_index(lower_bound);
#ifndef NDEBUG
		if (!__gij_of_r_numerator.count(atom_pair))
			throw Error("die : undefined __gij_of_r_numerator");
		if (idx >= __gij_of_r_numerator.at(atom_pair).size())
			throw Error("die : undefined __gij_of_r_numerator");
		if (idx >= __gij_of_r_bin_range_sum.size())
			throw Error("die : undefined __gij_of_r_bin_range_sum");
#endif
		double gij_of_r = __gij_of_r_numerator[atom_pair][idx] / __sum_gij_of_r_numerator[atom_pair];
		double denominator = __gij_of_r_bin_range_sum[idx] / __prot_lig_pairs.size();
		double ratio = gij_of_r / denominator;
		return -log(ratio);
	}
	double Score::__energy_cumulative(const pair_of_ints &atom_pair, const double &lower_bound) {
		const int idx = __get_index(lower_bound);
#ifndef NDEBUG
		if (!__gij_of_r_numerator.count(atom_pair))
			throw Error("die : undefined __gij_of_r_numerator");
		if (idx >= __gij_of_r_numerator.at(atom_pair).size())
			throw Error("die : undefined __gij_of_r_numerator");
		if (idx >= __bin_range_sum.size())
			throw Error("die : undefined __bin_range_sum");
#endif
		double numerator = __gij_of_r_numerator[atom_pair][idx] / __sum_gij_of_r_numerator[atom_pair];
		double denominator = __bin_range_sum[idx] / __total_quantity;
		double ratio = numerator / denominator;
		return -log(ratio);
	}
	double Score::non_bonded_energy(const Molecule &ligand) const {
		double energy_sum = 0.0;
		for (auto &assembly2 : ligand) {
		for (auto &model2 : assembly2) {
		for (auto &chain2 : model2) {
		for (auto &residue2 : chain2) {
		for (auto &atom2 : residue2) {
			const auto &atom_2 = atom2.idatm_type();
			dbgmsg("ligand atom = " << atom2.atom_number() << " crd= " << atom2.crd().pdb());
			for (auto &atom1 : __gridrec.get_neighbors(atom2, __dist_cutoff)) {
				const double dist = atom1->crd().distance(atom2.crd());
				const auto &atom_1 = atom1->idatm_type();
				const pair_of_ints atom_pair = minmax(atom_1, atom_2);
				const int idx = __get_index(dist);
				energy_sum += __energies.at(atom_pair).at(idx);
#ifndef NDEBUG				
				dbgmsg("ligand atom = " << atom2.atom_number() 
					<< " crd= " << atom2.crd().pdb() << "protein atom = " << atom1->atom_number() 
					<< " crd= " << atom1->crd().pdb() << " dist= " << dist 
					<< " atom_pair={" << atom_pair.first << "," << atom_pair.second << "}" 
					<< " lower_bound= " << __get_lower_bound(idx) << " energy_sum=" << energy_sum);
#endif
			}
		}}}}}
		return energy_sum;
	}
	double Score::non_bonded_energy(const AtomToCrd &atom_crd) const {
		double energy_sum = 0.0;
		for (auto &kv : atom_crd) {
			const Atom &atom2 = *kv.first;
			const Geom3D::Coordinate &atom2_crd = kv.second;
			const auto &atom_2 = atom2.idatm_type();
			dbgmsg("ligand atom = " << atom2.atom_number() << " crd= " << atom2_crd.pdb());
			for (auto &atom1 : __gridrec.get_neighbors(Atom(atom2_crd), __dist_cutoff)) {
				const double dist = atom1->crd().distance(atom2_crd);
				const auto &atom_1 = atom1->idatm_type();
				const pair_of_ints atom_pair = minmax(atom_1, atom_2);
				const int idx = __get_index(dist);
				energy_sum += __energies.at(atom_pair).at(idx);
#ifndef NDEBUG				
				dbgmsg("ligand atom = " << atom2.atom_number() 
					<< " crd= " << atom2_crd.pdb() << "protein atom = " << atom1->atom_number() 
					<< " crd= " << atom1->crd().pdb() << " dist= " << dist 
					<< " atom_pair={" << atom_pair.first << "," << atom_pair.second << "}" 
					<< " lower_bound= " << __get_lower_bound(idx) << " energy_sum=" << energy_sum);
#endif
			}
		}
		return energy_sum;
	}
	cluster::MapD<Molib::Molecule> Score::many_ligands_score(const Molib::Molecules &ligands) const {
		Benchmark::reset();
		cluster::MapD<Molib::Molecule> scores;
		cout << "Beginning scoring ..." << endl;
		for (auto &ligand : ligands) {
			dbgmsg("non_bonded_energy = " << non_bonded_energy(ligand)) ;
			scores[&ligand] = non_bonded_energy(ligand);
		}
		cout << "total time required for scoring " << ligands.size() 
			<< " decoys was " << Benchmark::seconds_from_start() << " wallclock seconds\n";
		return scores;
	}
	AtomTypeToEnergyPoint Score::compute_energy_grid(const set<int> &ligand_idatm_types, 
		const Geom3D::PointVec &gridpoints) {

		AtomTypeToEnergyPoint attep;
		Benchmark::reset();
		cout << "beginning scoring ..." << endl;
		for (const auto &ligand_atom_type : ligand_idatm_types) {
#ifndef NDEBUG
			int num = 0;
			dbgmsg("ligand_atom_type = " << help::idatm_unmask[ligand_atom_type]);
#endif
			for (auto &gridpoint : gridpoints) {
				double energy_sum = __distances_and_scores(gridpoint, ligand_atom_type);
				attep[ligand_atom_type].push_back({gridpoint, energy_sum});
				dbgmsg("compute energy grid : " << help::idatm_unmask[ligand_atom_type] 
					<< "\t" << gridpoint.with_underscores() << "\t" << setw(9) << setprecision(5) 
					<< fixed << energy_sum);
			}
		}
		cout << "total time required for scoring decoys was " 
			<< Benchmark::seconds_from_start() << " wallclock seconds\n";
		cout << "cleaning up temporary files ...\n";
		return attep;
	}
	AtomTypeToEnergyPoint Score::parse_energy_grid_file(const string &egrid_file) {
		AtomTypeToEnergyPoint attep;
		vector<string> scores;
		inout::Inout::read_file(egrid_file, scores);
		for (const string &line : scores) {
			boost::smatch m;
			if (boost::regex_search(line, m, boost::regex("^(\\S+)\\s+(\\S+)_(\\S+)_(\\S+)\\s+(\\S+)$"))) {
				if (m[1].matched && m[2].matched && m[3].matched && m[4].matched && m[5].matched) {
					attep[help::idatm_mask.at(m[1].str())].push_back(make_pair(Geom3D::Coordinate(stod(m[2].str()), 
						stod(m[3].str()), stod(m[4].str())), stod(m[5].str())));
					dbgmsg("processing atom type " << m[1].str());
				}
			}
		}
		return attep;
	}
};

ostream& operator<<(ostream& os, const Molib::AtomTypeToEnergyPoint &attep) {
	for (auto &kv : attep) {
		const int ligand_atom_type = kv.first;
		for (auto &ep : kv.second) {
			const Geom3D::Coordinate &point = ep.first;
			const double energy_sum = ep.second;
			os << help::idatm_unmask[ligand_atom_type] << "\t" 
				<< point.with_underscores() << "\t" << setw(9) 
				<< setprecision(5) << fixed << energy_sum << endl;
		}
	}
	return os;
}	

