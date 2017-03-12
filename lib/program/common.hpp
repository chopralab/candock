#ifndef COMMON_H
#define COMMON_H

#include "pdbreader/nrset.hpp"
#include "score/score.hpp"
#include "graph/graph.hpp"
#include "geom3d/geom3d.hpp"
#include "geom3d/coordinate.hpp"
#include "cluster/optics.hpp"
#include "ligands/genlig.hpp"
#include <thread>
#include <mutex>

namespace common {
	void change_residue_name(Molib::Molecule &ligand, const string &resn);
	void change_residue_name(Molib::Molecule &ligand, std::mutex &mtx, int &ligand_cnt);
	void create_mols_from_seeds(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
};

#endif
