#include "candock/linker/linker.hpp"
#include <queue>
#include "candock/cluster/greedy.hpp"
#include "statchem/geometry/geometry.hpp"
#include "statchem/graph/mcqd.hpp"
#include "statchem/helper/array2d.hpp"
#include "statchem/helper/benchmark.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/modeler/modeler.hpp"
#include "statchem/molib/bond.hpp"
#include "statchem/molib/nrset.hpp"
#include "statchem/score/score.hpp"

using namespace std;

namespace candock {
namespace linker {

Linker::Linker(OMMIface::Modeler& modeler, const molib::Molecule& receptor,
               const molib::Molecule& ligand, const molib::NRset& top_seeds,
               const molib::Atom::Grid& gridrec, const score::Score& score,
               const bool cuda, const bool iterative, const double dist_cutoff,
               const double spin_degrees, const double tol_seed_dist,
               const double lower_tol_seed_dist,
               const double upper_tol_seed_dist, const int max_possible_conf,
               const int link_iter, const double clash_coeff,
               const double docked_clus_rad, const double max_allow_energy,
               const int max_num_possibles, const int max_clique_size,
               const int max_iterations_final, const int max_iterations_pre,
               const std::string& platform, const std::string& precision,
               const std::string& accelerators) {
    if (iterative) {
        l = new IterativeLinker(
            modeler, receptor, ligand, top_seeds, gridrec, score, dist_cutoff,
            spin_degrees, tol_seed_dist, lower_tol_seed_dist,
            upper_tol_seed_dist, max_possible_conf, link_iter, clash_coeff,
            docked_clus_rad, max_allow_energy, max_num_possibles,
            max_clique_size, max_iterations_final, max_iterations_pre, platform,
            precision, accelerators);
    } else {
        l = new StaticLinker(
            modeler, receptor, ligand, top_seeds, gridrec, score, dist_cutoff,
            spin_degrees, tol_seed_dist, lower_tol_seed_dist,
            upper_tol_seed_dist, max_possible_conf, link_iter, clash_coeff,
            docked_clus_rad, max_allow_energy, max_num_possibles,
            max_clique_size, max_iterations_final, max_iterations_pre, platform,
            precision, accelerators);
    }
}

DockedConformation::Vec Linker::link() {
    l->init_openmm();
    Partial::Vec partial_conformations = l->init_conformations();
    return l->compute_conformations(partial_conformations);
}
};
}
