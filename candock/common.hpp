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
	class Centroid {
		Geom3D::Coordinate __centroid;
		double __radial_check;
	public:
		Centroid() {}
		Centroid(const Geom3D::Coordinate centroid, const double radial_check) :
			__centroid(centroid), __radial_check(radial_check) {}
		Geom3D::Coordinate get_centroid() const { return __centroid; }
		double get_radial_check() const { return __radial_check; }
	};

	class HCPoint {
		unique_ptr<Geom3D::Coordinate> __crd;
		double __energy;
		double __d;
	public:
		HCPoint(const Geom3D::Coordinate &c, const double &energy) : __energy(energy), 
			__crd(std::move(unique_ptr<Geom3D::Coordinate>(new Geom3D::Coordinate(c)))) {}
		Geom3D::Coordinate& crd() const { return *__crd; }
		double energy() const { return __energy; }
		double distance() const { return __d; } // NOT just dummy : needed by grid
		void distance(double d) { __d = d; } // NOT just dummy : needed by grid
		friend ostream& operator<< (ostream& stream, const HCPoint& h) {
			stream << setw(6) << left << "HETATM";
			stream << setw(5) << right << 1;
			stream << setw(1) << " ";
			stream << setw(4) << left << "X";
			stream << setw(1) << " ";
			stream << setw(3) << right << "DUM";
			stream << setw(1) << " ";
			stream << setw(1) << "A";
			stream << setw(4) << right << 1;
			stream << setw(1) << " ";
			stream << setw(27) << h.__crd->pdb();
			stream << setw(6) << setprecision(2) << fixed << right << 1.0;
			stream << setw(6) << setprecision(2) << fixed << right << h.__energy;
			stream << endl;
			return stream;
		}
	};
	
	typedef map<int, vector<unique_ptr<HCPoint>>> HCPoints;

	class PVertex : public template_vector_container<PVertex*, PVertex> {
		Molib::Atom &__v1;
		HCPoint &__v2;
		bool __visited;
		int __weight;
		Glib::node_id __index;
	public:
		PVertex(Molib::Atom &v1, HCPoint &v2, int ii) : __index(ii), __weight(int((10 - v2.energy())*10)), __v1(v1), __v2(v2), __visited(false) {}
		~PVertex() { dbgmsg("calling PVertex destructor"); }
		Geom3D::Coordinate& crd() const { return __v2.crd(); } // grid is done according to the second point (HCPoint)
		Molib::Atom &v1() const { return __v1; }
		const HCPoint &v2() const { return __v2; }
		double distance() const { return __v2.distance(); } // NOT just dummy : needed by grid
		void distance(double d) { __v2.distance(d); } // NOT just dummy : needed by grid
		void set_visited(bool b) { __visited = b; }
		bool visited() const { return __visited; }
		int weight() const { return __weight; }
		Glib::node_id get_index() const { return __index; }
		string get_label() const { return ""; } // dummy for graph ostream operator
		friend ostream& operator<< (ostream& stream, const PVertex& p) {
			return stream;
		}
	};

	typedef Glib::Graph<PVertex> ProductGraph;

	pair<cluster::Clusters<Molib::Molecule>, cluster::Clusters<Molib::Molecule>>
		cluster_molecules(const Molib::Molecules &mols,	const Molib::Score &score, 
		const double clus_rad, const int min_pts, const int max_num_clus, 
		const int max_mols_to_cluster=999999);

	vector<Centroid> set_centroids(const genlig::BindingSiteClusters &binding_site_clusters, 
		const double min_radial_check);
	vector<Centroid> set_centroids(const string &centroid_file, const double def_radial_check,
		const int num_bsites);
	Geom3D::PointVec identify_gridpoints(const Molib::Molecule &molecule, const Geom3D::Coordinate &centroid, 
		Molib::MolGrid &grid, const double &radial_check, const double &grid_spacing, const int &dist_cutoff,
		const double &excluded_radius, const double &max_interatomic_distance);
	HCPoints filter_scores(Molib::AtomTypeToEnergyPoint &attep, const double &top_percent);
	void create_mols_from_seeds(set<string> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
//~ #ifndef NDEBUG	
	void create_mols_from_fragments(set<string> &added, Molib::Molecules &seeds, const Molib::Molecules &mols);
//~ #endif
	ProductGraph product_graph(HCPoints &hcp, const Molib::Molecule &molecule, const double &grid_spacing);
	Molib::Molecules superimpose(ProductGraph::Cliques &maxclq, const Molib::Molecule &molecule);
	Molib::Molecules filter_clashes(const Molib::Molecules &rot_seeds, Grid<Molib::Atom> &gridrec);
	cluster::PairwiseDistances<Molib::Molecule> all_all_rmsd(const vector<Molib::Molecule*> &mols, const double &clus_rad);
	void convert_clusters_to_mols(Molib::Molecules &rep_mols, const cluster::Clusters<Molib::Molecule> &representatives);

	ostream& operator<<(ostream& os, const vector<Centroid>& centroids);
};

ostream& operator<<(ostream& os, const cluster::MapD<Molib::Molecule>& scores);
ostream& operator<<(ostream& os, const cluster::Clusters<Molib::Molecule>& molclus);
