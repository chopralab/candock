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
	Molib::NRset read_top_seeds_files(const Molib::Molecule &ligand, const string &top_seeds_dir, const string &top_seeds_file, const double top_percent);
	void create_mols_from_seeds(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
//~ #ifndef NDEBUG	
	void create_mols_from_fragments(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
//~ #endif
};
