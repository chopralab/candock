/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "linker.hpp"
#include "score/score.hpp"
#include "molib/nrset.hpp"
#include "molib/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "graph/mcqd.hpp"
#include "modeler/modeler.hpp"
#include "geom3d/geom3d.hpp"
#include "cluster/greedy.hpp"
#include <queue>

using namespace std;

namespace Linker {

	Linker::Linker(OMMIface::Modeler &modeler, const Molib::Molecule &receptor, 
			const Molib::Molecule &ligand, const Molib::NRset &top_seeds, 
			const Molib::Atom::Grid &gridrec, const Molib::Score &score,
			const bool cuda, const bool iterative, const double dist_cutoff, 
			const double spin_degrees, const double tol_seed_dist, 
			const double lower_tol_seed_dist, const double upper_tol_seed_dist, 
			const int max_possible_conf, const int link_iter,
			const double clash_coeff, const double docked_clus_rad,
			const double max_allow_energy, const int max_num_possibles, 
			const int max_clique_size, const int max_iterations_final) { 

		
		if (iterative) {
			l = new IterativeLinker(modeler, receptor, ligand, top_seeds, 
				gridrec, score, dist_cutoff, spin_degrees, tol_seed_dist, 
				lower_tol_seed_dist, upper_tol_seed_dist, max_possible_conf, link_iter, 
				clash_coeff, docked_clus_rad, max_allow_energy, max_num_possibles, 
				max_clique_size, max_iterations_final);
		} else {
			l = new StaticLinker(modeler, receptor, ligand, top_seeds, 
				gridrec, score, dist_cutoff, spin_degrees, tol_seed_dist, 
				lower_tol_seed_dist, upper_tol_seed_dist, max_possible_conf, link_iter, 
				clash_coeff, docked_clus_rad, max_allow_energy, max_num_possibles, 
				max_clique_size, max_iterations_final);
		}
	}
	
	DockedConformation::Vec Linker::link() {
		l->init_openmm();
		Partial::Vec partial_conformations = l->init_conformations();
		return l->compute_conformations(partial_conformations);
	}

};
