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
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

ostream& operator<<(ostream& os, const cluster::PairwiseDistances<Molib::Molecule>& pairwise_distances)	{
	for (auto &kv1 : pairwise_distances) {
		const Molib::Molecule &molecule1 = *kv1.first;
		for (auto &kv2 : kv1.second) {
			const Molib::Molecule &molecule2 = *kv2.first;
			const double &distance = kv2.second;
			os << molecule1.name() << "\t" << molecule1.first().first().number() 
				<< "\t" << molecule2.first().first().number() 
				<< "\t" << distance << endl;
		}
	}
	return os;
}	

ostream& operator<<(ostream& os, const cluster::MapD<Molib::Molecule>& scores)	{
	for (auto &kv : scores) {
		const Molib::Molecule &molecule = *kv.first;
		const double &energy = kv.second;
		os << molecule.name() << "\t" << molecule.first().first().number() 
			<< "\t" << setprecision(5) << fixed << energy << endl;
	}
	return os;
}

ostream& operator<<(ostream& os, const cluster::Clusters<Molib::Molecule>& molclus)	{
	for (auto &kv : molclus) {
		const int &cluster_number = kv.first;
		Molib::Molecule &member = *kv.second;
		member.set_name(member.name() + "_" + help::to_string(cluster_number));
		os << member;
	}
	return os;
}	

namespace common {



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

	/* Part7a stuff
	 * 
	 */
	Molib::NRset read_top_seeds_files(const Molib::Molecule &ligand, const string &top_seeds_file) {
		Molib::NRset top_seeds;
		const Molib::Model &model = ligand.first().first();
		for (auto &kv : model.get_seeds()) { // iterate over seeds
			const string &nm = kv.first;
			dbgmsg("reading top_seeds_file for seed number = " << nm);
			Molib::PDBreader pdb("tmp/" + nm + "/" + top_seeds_file, 
				Molib::PDBreader::all_models);
			top_seeds.add(new Molib::Molecules(pdb.parse_molecule()));
		}
		return top_seeds;
	}
	void create_mols_from_seeds(set<string> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
		for (auto &molecule : mols)
		for (auto &assembly : molecule)
		for (auto &model : assembly) {
			for (auto &kv : model.get_seeds()) { // iterate over seeds
				const string &nm = kv.first;
				dbgmsg("considering to add " << nm);
				if (!added.count(nm)) { // take seeds that haven't been docked already
					dbgmsg("added " << nm);
					added.insert(nm);
					if (kv.second.empty()) throw Error("die : no seed exists with this name : " + nm);
					const Molib::AtomSet &seed_atoms = *kv.second.begin();
					// add to new molecules
					Molib::Molecule &seed = seeds.add(new Molib::Molecule(nm));
					Molib::Assembly &a = seed.add(new Molib::Assembly(0));
					Molib::Model &mod = a.add(new Molib::Model(1));
					Molib::Chain &c = mod.add(new Molib::Chain('X'));
					Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
					for (const Molib::Atom *atom : seed_atoms) {
						Molib::Atom &at = r.add(new Molib::Atom(*atom));
						dbgmsg("added atom = " << at);
					}
					seed.regenerate_bonds(molecule);
				}
			}
		}
	}
//~ #ifndef NDEBUG
	void create_mols_from_fragments(set<string> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
		for (auto &molecule : mols)
		for (auto &assembly : molecule)
		for (auto &model : assembly) {
			for (auto &kv : model.get_rigid()) { // iterate over rigid fragments
				const string &nm = kv.first;
				dbgmsg("considering to add " << nm);
				if (!added.count(nm)) { // take seeds that haven't been docked already
					dbgmsg("added " << nm);
					added.insert(nm);
					if (kv.second.empty()) throw Error("die : no seed exists with this name : " + nm);
					const Molib::AtomSet &seed_atoms = *kv.second.begin();
					// add to new molecules
					Molib::Molecule &seed = seeds.add(new Molib::Molecule(nm));
					Molib::Assembly &a = seed.add(new Molib::Assembly(0));
					Molib::Model &mod = a.add(new Molib::Model(1));
					Molib::Chain &c = mod.add(new Molib::Chain('X'));
					Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
					for (const Molib::Atom *atom : seed_atoms) {
						Molib::Atom &at = r.add(new Molib::Atom(*atom));
						dbgmsg("added atom = " << at);
					}
					seed.regenerate_bonds(molecule);
				}
			}
		}
	}
//~ #endif

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
			vector<Molib::MolGraph> mol_graphs;
			/* compute geometric centers of mols
			 */
			for (int i = 0; i < mols.size(); ++i) {
				int sz = 0;
				mol_graphs.push_back(create_graph(mols[i]->get_atoms()));
				for (auto &patom : mols[i]->get_atoms()) {
					geom_centers[i] = geom_centers[i] + patom->crd();
					++sz;
				}
				geom_centers[i] = geom_centers[i] / sz;
			}
			/* compute matches between same graphs
			 */
			Molib::MolGraph::Matches matches = mol_graphs[0].match(mol_graphs[0]);
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
};
