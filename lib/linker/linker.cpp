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
