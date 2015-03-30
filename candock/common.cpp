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

	//~ ostream& operator<<(ostream& os, const vector<common::Centroid>& centroids) {
	ostream& operator<<(ostream& os, const vector<Centroid>& centroids) {
		for (auto &centroid : centroids) {
			os << centroid.get_centroid().simple() << " " << fixed
				<< setprecision(3) << centroid.get_radial_check() << endl;
		}
		return os;
	}
	


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
		//~ cluster_molecules(const Molib::Molecules &mols, Molib::Molecules &rep_mols,
		cluster_molecules(const Molib::Molecules &mols,	const Molib::Score &score, 
		const double clus_rad, const int min_pts, const int max_num_clus, const int max_mols_to_cluster) {

		// if only one molecule, just return it ... (cluster algorithm
		// requires at least two molecules due to pairwise_distances)
		//~ if (mols.size() == 1)
			//~ return cluster::Clusters<Molib::Molecule> {{1, &mols.first()}};
		if (mols.size() == 1)
			return make_pair(cluster::Clusters<Molib::Molecule> {{1, &mols.first()}},
				cluster::Clusters<Molib::Molecule> {{1, &mols.first()}});
		cluster::MapD<Molib::Molecule> scores = 
			score.many_ligands_score(mols);
//~ //////////////////////
		//~ cluster::MapD<Molib::Molecule> scores;
		//~ double sc = -100;
		//~ for (auto &molecule : mols) {
			//~ scores.insert({&molecule, sc});
			//~ sc += 200;
		//~ }
//~ ////////////////////////
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
		//~ return optics.extract_dbscan(clus_rad, max_num_clus, true).second;
		return optics.extract_dbscan(clus_rad, max_num_clus, true);
		//~ auto clusters_reps_pair = optics.extract_dbscan(clus_rad, max_num_clus, true);
		//~ convert_clusters_to_mols(rep_mols, clusters_reps_pair.second);
		//~ if (1 == 1) {
			//~ inout::output_file(clusters_reps_pair.first, "clusters.txt", ios_base::app); // get the PDB file with all the docked fragment clusters
			//~ throw Error("finished : getting docked fragment clusters...exiting...");
		//~ }
	}

	/* Centroid stuff
	 * 
	 */
	vector<Centroid> split_binding_site(const Molib::Molecules &binding_site_ligands, 
		const double def_radial_check) {

		vector<Centroid> centroids;
		Molib::AtomVec atoms;
		for (auto &molecule : binding_site_ligands) {
			auto molecule_atoms = molecule.get_atoms();
			atoms.insert(atoms.end(), molecule_atoms.begin(), molecule_atoms.end());
		}
		gsl_matrix *data = gsl_matrix_alloc(3, atoms.size());
		for (int i = 0; i < atoms.size(); ++i) {
			gsl_matrix_set(data, 0, i, atoms[i]->crd().x());
			gsl_matrix_set(data, 1, i, atoms[i]->crd().y());
			gsl_matrix_set(data, 2, i, atoms[i]->crd().z());
		}
		gsl_matrix *projection, *eigenmatrix;
		gsl_vector *mean;
		// get the first principal component only
		tie(projection, eigenmatrix, mean) = Geom3D::pca(data, 1);
		const Geom3D::Point geom_center = *mean;
		const Geom3D::Vector3 unit_vector = gsl_matrix_column(eigenmatrix, 0).vector;
		vector<double> projected_points;
		const unsigned int cols = projection->size2;
		for (int i = 0; i < cols; ++i) 
			projected_points.push_back(gsl_matrix_get(projection, 0, i)); 
		auto ret = minmax_element(projected_points.begin(), projected_points.end());
		const double min_value = *ret.first, max_value = *ret.second;
		const double max_dist = max_value - min_value;
		const int num_centroids = ceil(max_dist / def_radial_check);
		//~ const int num_centroids = 1 + ceil(max_dist / def_radial_check); // add one just in case...
		dbgmsg("num_centroids = " << num_centroids << " min_value = " << min_value
			<< " max_dist = " << max_dist);
		const double interval = max_dist / (num_centroids + 1);
		for (int i = 1; i <= num_centroids; ++i) {
			const double position = min_value + i * interval;
			centroids.push_back(Centroid(Geom3D::line_evaluate(
				geom_center, unit_vector, position), def_radial_check));
			dbgmsg("original centroid geom_center = " << geom_center
				<< " unit_vector = " << unit_vector << " position = "
				<< position << " : adding splitted centroid[" << i << "] at center = "
				<< Geom3D::line_evaluate(geom_center, unit_vector, position));
			dbgmsg("ATOM      1  X   GEO X   1    " 
				<< Geom3D::line_evaluate(geom_center, unit_vector, position).pdb()
				<< "  1.00  0.00           N  ");
		}
		gsl_matrix_free(data);
		gsl_matrix_free(projection);
		gsl_matrix_free(eigenmatrix);
		gsl_vector_free(mean);
		return centroids;
	}

	vector<Centroid> set_centroids(const genlig::BindingSiteClusters &binding_site_clusters, 
		const double def_radial_check) {
		vector<Centroid> centroids;
		if (binding_site_clusters.empty()) 
			throw Error("die : no binding sites could be predicted for this protein \
(set binding site(s) centroids manually through centroid file)");
		for (auto &kv : binding_site_clusters) {
			const int cluster_number = kv.first;
			const Molib::Molecules &binding_site_ligands = kv.second;
			double radial_check = 2 * binding_site_ligands.compute_max_radius();
			centroids.push_back(Centroid(binding_site_ligands.compute_geometric_center(), 
				radial_check));
			//~ if (radial_check <= def_radial_check) {
				//~ cerr << "note : radial check calculated using ProBiS " 
					//~ << " binding site detection is too small (" << radial_check
					//~ << ") using default value of " << def_radial_check << endl;
				//~ radial_check = def_radial_check;
				//~ centroids.push_back(Centroid(binding_site_ligands.compute_geometric_center(), 
					//~ radial_check));
			//~ } else {
				//~ auto result = split_binding_site(binding_site_ligands, def_radial_check);
				//~ cerr << "note : radial check calculated using ProBiS " 
					//~ << " binding site detection is too LARGE (" << radial_check
					//~ << ") splitting to " << result.size() << " smaller binding sites of "
					//~ << def_radial_check << " Angstroms." << endl;
				//~ centroids.insert(centroids.end(), result.begin(), result.end()); 
			//~ }
			dbgmsg("setting centroids from predicted binding sites"
				<< " cluster number " << cluster_number 
				<< " radial check = " << centroids.back().get_radial_check() 
				<< " centroid = " << centroids.back().get_centroid());
		}
		dbgmsg(centroids);
		return centroids;
	}
	vector<Centroid> set_centroids(const string &centroid_file, const double def_radial_check,
		const int num_bsites) {
			
		vector<Centroid> centroids;
		vector<string> data;
		inout::Inout::read_file(centroid_file, data);
		for (string &line : data) {
			stringstream ss(line);
			double x, y, z, rc = -1;
			ss >> x >> y >> z;
			ss >> rc;
			if (rc == -1) {
				rc = def_radial_check;
				cerr << "warning : radial check is missing in the centroid file - "
					<< "using default value of " << rc << endl;
			}
			centroids.push_back(Centroid(Geom3D::Coordinate(x, y, z), rc));
			dbgmsg("setting centroids from file"
				<< " number " << centroids.size() 
				<< " radial check = " << centroids.back().get_radial_check() 
				<< " centroid = " << centroids.back().get_centroid());
			if (centroids.size() >= num_bsites)
				break;
		}
		if (centroids.empty()) 
			throw Error("die: could not find centroid in centroid file " 
				+ centroid_file + "\n");
		return centroids;
	}

	/* Part1 stuff
	 * 
	 */
	Geom3D::PointVec identify_gridpoints(const Molib::Molecule &molecule, const Geom3D::Coordinate &centroid, 
		Molib::MolGrid &grid, const double &radial_check, const double &grid_spacing, const int &dist_cutoff,
		const double &excluded_radius, const double &max_interatomic_distance) {
			
		Geom3D::PointVec gridpoints;
		Benchmark::reset();
		Geom3D::Coordinate min = centroid - ceil(radial_check);
		Geom3D::Coordinate max = centroid + ceil(radial_check);
		const int total_gridpoints = 3*ceil((max.x()-min.x())/grid_spacing)
									*ceil((max.y()-min.y())/grid_spacing)
									*ceil((max.z()-min.z())/grid_spacing);
		cout <<  "approximately " << total_gridpoints << " gridpoints to evaluate\n\n";
		int model_number = 1;
		int points_kept = 0;
		int gridpoint_counter = 0;
		const double r = grid_spacing/2;
		const int last_column = ceil(2*radial_check/r);
		const int last_row = ceil(2*radial_check/(sqrt(3)*r));
		const int last_layer = ceil(2*radial_check/(2*r*sqrt(6)/3));
		Geom3D::Coordinate eval;
		for(int column=0;column<=last_column;column++) {
			int even_column=(column%2==0) ? 1 : 0; // 1 if odd, 0 if even
			for(int row=0;row<=last_row;row++) {
				int even_row=(row%2==0) ? 1 : 0; // 1 if odd, 0 if even
				for(int layer=0;layer<=last_layer;layer++) {
					int even_layer=(layer%2==0) ? 1 : 0; // 1 if odd, 0 if even	
					if ((even_column==0 && even_row==0) || (even_column==1 && even_row==1)) {
						if (even_layer==1) {
							eval.set_x(min.x()+column*r);
							eval.set_y(min.y()+sqrt(3)*r*row);
							eval.set_z(min.z()+layer*2*r*sqrt(6)/3);
						}
						else {
							eval.set_x(min.x()+r+column*r);		 
							eval.set_y(min.y()+r/sqrt(3)+sqrt(3)*r*row); 
							eval.set_z(min.z()+layer*2*r*sqrt(6)/3);      
						}
						int okay_min=1;
						int okay_max=1;
						double closest = 10000.0;
						const double radial_dist = eval.distance(centroid);
						if (radial_dist <= radial_check) {
							Molib::Atom at(eval);
							vector<Molib::Atom*> neighbors = grid.get_neighbors(at, dist_cutoff);
							for (Molib::Atom *a : neighbors) {
								Molib::Atom &atom = *a;
								const double vdW = help::vdw_radius[atom.idatm_type()];
								dbgmsg("vdW = " << vdW);
								const double eval_dist = excluded_radius+0.9*vdW;
								const double distance = atom.crd().distance(eval);
								if (distance <= eval_dist) {
									okay_min = 0;
									dbgmsg("distance = " << distance << " eval_dist = " << eval_dist);
									goto OUTER;
								}
								else {
									okay_min = 1;
									if (distance < closest) closest = distance;
								}
							}
							OUTER:
							if (closest>max_interatomic_distance) okay_max = 0;
							if (okay_min*okay_max > 0) {
								gridpoints.push_back(eval);
								model_number++;
								points_kept++;
							}
						}
					}
					const int mod = gridpoint_counter % 10000;
					if (mod == 0) {
						double wall_secs = Benchmark::seconds_from_start();
						cout << "Processing gridpoint " << gridpoint_counter 
							<< " of approximately " << total_gridpoints 
							<< " (took " << wall_secs << " seconds)\n";
						Benchmark::reset();
					}
					gridpoint_counter++;
				}
				dbgmsg("column = " << column);
			}
		}
		// the last ones that did not get to the next mod==0
		cout << points_kept << " points kept out of " << gridpoint_counter 
			<< " total gridpoints\n";
		return gridpoints;
	}
	
	/* Part7a stuff
	 * 
	 */
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
	HCPoints filter_scores(Molib::AtomTypeToEnergyPoint &attep, const double &top_percent) {
		HCPoints hcp;
		for (auto &u : attep) {
			sort(u.second.begin(), u.second.end(), [] (const pair<Geom3D::Coordinate, double> &i, 
				const pair<Geom3D::Coordinate, double> &j) { return i.second < j.second; });
		}
#ifndef NDEBUG
		for (auto &u : attep) for (auto &i : u.second) { dbgmsg(u.first << " " << i.first << " " << i.second); }
#endif
		for (auto &u : attep) {
			const vector<pair<Geom3D::Coordinate, double>> &points = u.second;
			const int number = ceil(top_percent * points.size()); // return top_percent, e.g. 40% of top scores
			if (points.size() < number) 
				throw Error("die : Invalid top_percent given (should be between 0 and 1)");
			const auto &atom_type = u.first;
			for (int i = 0; i < number; i++) {
				const Geom3D::Coordinate &crd = points[i].first;
				const double energy = points[i].second;
				hcp[atom_type].push_back(unique_ptr<HCPoint>(new HCPoint(crd, energy)));
			}
		}
#ifndef NDEBUG
		for (auto &kv1 : hcp) for (auto &kv2 : kv1.second) { dbgmsg("atom type = " << kv1.first << " " << *kv2); }
		for (auto &kv1 : hcp) {
			stringstream ss;
			for (auto &kv2 : kv1.second) {
				ss << *kv2;
			}
			inout::Inout::file_open_put_stream("topscores_" + 
				help::to_string(help::idatm_unmask[kv1.first]) + ".pdb", ss);
		}
#endif
		return hcp;
	}
	ProductGraph product_graph(HCPoints &hcp, const Molib::Molecule &molecule, const double &grid_spacing) {
		Benchmark::reset();
		cout << "set up product graph" << endl;
		dbgmsg("for molecule = " << endl << molecule);
		int i = 0;
		vector<unique_ptr<PVertex>> vertices;
		for (auto &patom : molecule.get_atoms()) {
			auto &atom = *patom;
			dbgmsg(atom);
			dbgmsg(atom.idatm_type() << " " << (!hcp.count(atom.idatm_type()) ? "error" : "cool"));
			for (auto &point : hcp[atom.idatm_type()]) {
				vertices.push_back(unique_ptr<PVertex>(new PVertex(atom, *point, i++)));
			}
		}
		dbgmsg("out of set up product graph");
		ProductGraph graph(std::move(vertices), false, false); // adjacency matrix will come later...
		dbgmsg("number of vertices = " << i);
		cout << "time to generate vertices of product graph " << Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		vector<PVertex*> pvertices;
		for (auto &pv : graph) {
			pvertices.push_back(&pv);
		}
		Grid<PVertex> pgrid(pvertices);
		cout << "time to make grid from product graph " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		// calculate max distance between atoms of seed
		const double neighb_dist = molecule.max_dist() + grid_spacing;
		dbgmsg("taking neighbor distance " << neighb_dist);
		graph.init_conn(); // resize product graph conn table
		for (auto &pv1 : graph) {
			vector<PVertex*> neighbors = pgrid.get_neighbors(pv1, neighb_dist);
			pv1.set_visited(true);
			for (auto &pv2 : neighbors) {
				if (!pv2->visited() && &pv1.v1() != &pv2->v1() && &pv1.v2() != &pv2->v2()) {
					const double dist1 = pv1.v1().crd().distance_sq(pv2->v1().crd());
					const double dist2 = pv2->v2().distance();
					const double term = dist1 + dist2 - grid_spacing * grid_spacing;
					if ( term * term < 4 * dist1 * dist2) {
						pv1.add(pv2);
						pv2->add(&pv1);
						graph.set_conn(pv1.get_index(), pv2->get_index());
					}
				}
			}
		}
#ifndef NDEBUG
		stringstream ss;
		ss << graph;
		inout::Inout::file_open_put_stream("pgraph.dbg", ss);
#endif
		cout << "time to generate edges of product graph " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		return graph;
	}
	Molib::Molecules superimpose(ProductGraph::Cliques &maxclq, const Molib::Molecule &molecule) {
		Benchmark::reset();
#ifndef NDEBUG
		for (auto &clique: maxclq) { for (auto &vertex : clique) {	dbgmsg(((PVertex*)vertex)->v1()); } dbgmsg("TER"); }
		for (auto &clique: maxclq) { for (auto &vertex : clique) {	dbgmsg(((PVertex*)vertex)->v2()); } dbgmsg("TER"); }
#endif
		Molib::Molecules rot_mols;
		for (auto &clique : maxclq) {
			try {
				dbgmsg("clique.size = " << clique.size());
				if (clique.size() > 2) {
					Kabsch k;
					k.resize(clique.size());
					for (auto &vertex : clique) {
						// align left coords to right coords
						k.add_vertex(((PVertex*)vertex)->v1().crd(), 
							((PVertex*)vertex)->v2().crd());
					}
					dbgmsg("before superimpose");
					k.superimpose();
					const Geom3D::Matrix rm = k.get_rota();
					dbgmsg(rm);
					Molib::Molecule &added_mol = rot_mols.add(new Molib::Molecule(molecule)); // copy fragment to new molecule
					added_mol.rotate(rm); // rotate it
					dbgmsg(added_mol);
				} 
#ifndef NDEBUG
				else {
					dbgmsg("skipping superimposition due to clique size < 3 : clique size = " 
						<< clique.size());
				}
#endif
			} catch (exception& e) {
				cerr << "superimposition failed (nothing to worry about) clique size = " 
					<< clique.size() << " msg = " << e.what() << " ... skipping ..." << endl;
			}
		}
		cout << "time to superimpose " << Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		return rot_mols;
	}
	Molib::Molecules filter_clashes(const Molib::Molecules &rot_seeds, Grid<Molib::Atom> &gridrec) {
		Benchmark::reset();
		Molib::Molecules non_clashing_seeds;
		int i=0;
		for (auto &molecule : rot_seeds) {
			for (auto &patom : molecule.get_atoms()) 
				if(gridrec.has_neighbor_within(*patom, 1.0)) goto has_clashes;
			non_clashing_seeds.add(new Molib::Molecule(molecule));
			non_clashing_seeds.last().last().last().set_number(++i); // set number of model to the conformation of docked seed
			;
			has_clashes:
			;
		}
		cout << "time to filter clashes took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
		return non_clashing_seeds;
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
			vector<Molib::MolGraph> mol_graphs;
			/* compute geometric centers of mols
			 */
			for (int i = 0; i < mols.size(); ++i) {
				int sz = 0;
				//~ mol_graphs.push_back(create_graph(*mols[i]));
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
	//~ /* Part9 stuff
	 //~ * 
	 //~ */
	//~ void convert_clusters_to_mols(Molib::Molecules &rep_mols, 
		//~ const cluster::Clusters<Molib::Molecule> &representatives) {
//~ 
		//~ for (auto &kv : representatives) {			
			//~ Molib::Molecule &rep = *kv.second;
			//~ dbgmsg("representative of cluster number " << kv.first << " is molecule " 
				//~ << rep.name() << " conformation " << rep.first().first().number());
			//~ rep_mols.add(new Molib::Molecule(rep));
		//~ }
	//~ }
};
