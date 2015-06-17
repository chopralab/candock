#include "pdbreader/molecule.hpp"
#include "score/score.hpp"
#include "graph/graph.hpp"
#include "geom3d/geom3d.hpp"
#include "geom3d/coordinate.hpp"
#include "cluster/optics.hpp"
#include "ligands/genlig.hpp"
#include <thread>
#include <mutex>

namespace common {
	pair<cluster::Clusters<Molib::Molecule>, cluster::Clusters<Molib::Molecule>>
		cluster_molecules(const Molib::Molecules &mols,	const Molib::Score &score, 
		const double clus_rad, const int min_pts, const int max_num_clus, 
		const int max_mols_to_cluster=999999);

	Molib::NRset read_top_seeds_files(const Molib::Molecule &ligand, const string &top_seeds_file);
	void create_mols_from_seeds(set<string> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
//~ #ifndef NDEBUG	
	void create_mols_from_fragments(set<string> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
//~ #endif
	Molib::Molecules dock_seeds(Geom3D::GridPoints &gridpoints, const Molib::Molecule &molecule, const double &grid_spacing);
	cluster::PairwiseDistances<Molib::Molecule> all_all_rmsd(const vector<Molib::Molecule*> &mols, const double &clus_rad);
	void convert_clusters_to_mols(Molib::Molecules &rep_mols, const cluster::Clusters<Molib::Molecule> &representatives);

};

ostream& operator<<(ostream& os, const cluster::MapD<Molib::Molecule>& scores);
ostream& operator<<(ostream& os, const cluster::Clusters<Molib::Molecule>& molclus);
