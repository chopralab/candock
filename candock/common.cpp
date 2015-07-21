#include "common.hpp"
#include "helper/inout.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "geom3d/matrix.hpp"
#include "geom3d/pca.hpp"
#include "geom3d/geom3d.hpp"
#include "kabsch/kabsch.hpp"
#include "score/score.hpp"
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/atom.hpp"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

namespace common {

	/* Part7a stuff
	 * 
	 */
	Molib::NRset read_top_seeds_files(const Molib::Molecule &ligand, const string &top_seeds_file) {
		Molib::NRset top_seeds;
		const Molib::Model &model = ligand.first().first();
		for (auto &fragment : model.get_rigid()) { // iterate over seeds
			if (fragment.is_seed()) {
				dbgmsg("reading top_seeds_file for seed id = " << fragment.get_seed_id());
				Molib::PDBreader pdb("tmp/" + help::to_string(fragment.get_seed_id()) + "/" + top_seeds_file, 
					Molib::PDBreader::all_models);
				top_seeds.add(new Molib::Molecules(pdb.parse_molecule()));
			}
		}
		return top_seeds;
	}

	void create_mols_from_seeds(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
		for (auto &molecule : mols)
		for (auto &assembly : molecule)
		for (auto &model : assembly) {
			for (auto &fragment : model.get_rigid()) { // iterate over seeds
				if (fragment.is_seed()) {
					dbgmsg("considering to add " << fragment.get_seed_id());
					if (!added.count(fragment.get_seed_id())) { // take seeds that haven't been docked already
						dbgmsg("added " << fragment.get_seed_id());
						added.insert(fragment.get_seed_id());
						// add to new molecules
						Molib::Molecule &seed = seeds.add(new Molib::Molecule(help::to_string(fragment.get_seed_id())));
						Molib::Assembly &a = seed.add(new Molib::Assembly(0));
						Molib::Model &mod = a.add(new Molib::Model(1));
						Molib::Chain &c = mod.add(new Molib::Chain('X'));
						Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
						for (const Molib::Atom *atom : fragment.get_all()) {
							Molib::Atom &at = r.add(new Molib::Atom(*atom));
							dbgmsg("added atom = " << at);
						}
						seed.regenerate_bonds(molecule);
					}
				}
			}
		}
	}

//~ #ifndef NDEBUG
	void create_mols_from_fragments(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
		for (auto &molecule : mols)
		for (auto &assembly : molecule)
		for (auto &model : assembly) {
			for (auto &fragment : model.get_rigid()) { // iterate over seeds
				if (fragment.is_seed()) {
					dbgmsg("considering to add " << fragment.get_seed_id());
					if (!added.count(fragment.get_seed_id())) { // take seeds that haven't been docked already
						dbgmsg("added " << fragment.get_seed_id());
						added.insert(fragment.get_seed_id());
						// add to new molecules
						Molib::Molecule &seed = seeds.add(new Molib::Molecule(help::to_string(fragment.get_seed_id())));
						Molib::Assembly &a = seed.add(new Molib::Assembly(0));
						Molib::Model &mod = a.add(new Molib::Model(1));
						Molib::Chain &c = mod.add(new Molib::Chain('X'));
						Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
						for (const Molib::Atom *atom : fragment.get_all()) {
							Molib::Atom &at = r.add(new Molib::Atom(*atom));
							dbgmsg("added atom = " << at);
						}
						seed.regenerate_bonds(molecule);
					}
				}
			}
		}
	}
//~ #endif

	/* Cluster stuff
	 * 
	 */

	void convert_clusters_to_mols(Molib::Molecules &rep_mols, 
		const cluster::Clusters<Molib::Molecule> &representatives) {

		for (auto &kv : representatives) {			
			Molib::Molecule &rep = *kv.second;
			dbgmsg("representative of cluster number " << kv.first << " is molecule " 
				<< rep.name() << " conformation " << rep.first().first().number());
			rep_mols.add(new Molib::Molecule(rep));
		}
	}
	
	pair<cluster::Clusters<Molib::Molecule>, cluster::Clusters<Molib::Molecule>>
		cluster_molecules(const Molib::Molecules &mols,	const Molib::Score &score, 
		const double clus_rad, const int min_pts, const int max_num_clus, const int max_mols_to_cluster) {
		
		// if only one molecule, just return it ... (cluster algorithm
		// requires at least two molecules due to pairwise_distances)
		if (mols.size() == 1)
			return make_pair(cluster::Clusters<Molib::Molecule> {{1, &mols.first()}},
				cluster::Clusters<Molib::Molecule> {{1, &mols.first()}});
		
		cluster::MapD<Molib::Molecule> scores = 
			score.many_ligands_score(mols);
		
		dbgmsg(scores);
		vector<Molib::Molecule*> mols_to_cluster;
		set<double> top_scores;
		for (auto &kv : scores) top_scores.insert(kv.second);
		for (auto &kv : scores) {
			Molib::Molecule &molecule = *kv.first;
			const double score = kv.second;
			if (distance(top_scores.begin(), top_scores.find(score)) < max_mols_to_cluster) {
				mols_to_cluster.push_back(&molecule);
			}
		}
		
		cluster::PairwiseDistances<Molib::Molecule> pairwise_distances = 
			common::all_all_rmsd(mols_to_cluster, clus_rad);
		
		dbgmsg(pairwise_distances);
		cluster::Optics<Molib::Molecule, std::less<double>, std::greater<double>> optics(
			pairwise_distances, scores, clus_rad + 0.1, min_pts);
		
		return optics.extract_dbscan(clus_rad, max_num_clus, true);
	}

	/* Part8 stuff
	 * 
	 */
	cluster::PairwiseDistances<Molib::Molecule> all_all_rmsd(const vector<Molib::Molecule*> &mols, 
		const double &clus_rad) {

		Benchmark::reset();
		cluster::PairwiseDistances<Molib::Molecule>  pairwise_distances;
		cout << "calculating all-against-all rmsd" << endl;
		if (!mols.empty()) {
			vector<Geom3D::Coordinate> geom_centers(mols.size());
			vector<Molib::Atom::Graph> mol_graphs;
			/* compute geometric centers of mols
			 */
			for (int i = 0; i < mols.size(); ++i) {
				int sz = 0;
				mol_graphs.push_back(Molib::Atom::create_graph(mols[i]->get_atoms()));
				for (auto &patom : mols[i]->get_atoms()) {
					geom_centers[i] = geom_centers[i] + patom->crd();
					++sz;
				}
				geom_centers[i] = geom_centers[i] / sz;
			}
			/* compute matches between same graphs
			 */
			Molib::Atom::Graph::Matches matches = mol_graphs[0].match(mol_graphs[0]);
			const double clus_rad_sq = clus_rad * clus_rad;
			for (int i = 0; i < mols.size(); ++i) {
				for(int j = i + 1; j < mols.size(); ++j) {
					const double geom_dist_sq = 
						geom_centers[i].distance_sq(geom_centers[j]); // faster than rmsd
					const double rmsd = 
						(geom_dist_sq < clus_rad_sq ? 
						compute_rmsd(mol_graphs[i], mol_graphs[j], matches) 
						: sqrt(geom_dist_sq));
					pairwise_distances[mols[i]][mols[j]] = rmsd;
				}
			}
		}
		cout << "total time required for rmsd calculation was " 
			<< Benchmark::seconds_from_start() << " wallclock seconds\n";
		return pairwise_distances;
	}

	double compute_rmsd(const Molib::Atom::Graph &g1, const Molib::Atom::Graph &g2, 
		const Molib::Atom::Graph::Matches &m) {
		dbgmsg("calculate rmsd between two conformations of the same \
			molecule (can do symmetric molecules such as benzene, etc.)");
		double min_sum_squared = HUGE_VAL;
		// try calculating rmsd of each mapping of molecule to molecule ...
		for (auto &mv : m) {
			auto &vertices1 = mv.first;
			auto &vertices2 = mv.second;
			// match must stretch over the whole graph g1 (and g2)
			if (vertices2.size() == g1.size()) {
				double sum_squared = 0;
				for (int i = 0; i < vertices2.size(); ++i) {
					Molib::Atom &a1 = g1[vertices1[i]];
					Molib::Atom &a2 = g2[vertices2[i]];
					sum_squared += a1.crd().distance(a2.crd());
				}
				if (sum_squared < min_sum_squared)
					min_sum_squared = sum_squared;
			} else { 
				throw Error("die : RMSD calculation aborted due to imperfect match"); 
			}
		}
		if (min_sum_squared == HUGE_VAL) throw Error("die : RMSD could not be calculated");
		return sqrt(min_sum_squared / g1.size());	
	}
};
