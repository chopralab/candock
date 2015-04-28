#include <iostream>
#include <exception>
#include <typeinfo>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"
#include "geom3d/matrix.hpp"
#include "helper/error.hpp"
#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "cluster/optics.hpp"
#include "jsonreader.hpp"
#include "nosqlreader.hpp"
#include "common.hpp"
#include "genclus.hpp"
using namespace std;

namespace genclus {
	string remove_special_characters(string str) {
		str.erase(std::remove_if(str.begin(), str.end(), [](const char &c ) { return c==':'||c=='"'||c=='\\';}), str.end());
		return str;
	}
	void add_to_json(JsonReader &jr, cluster::Clusters<Molib::Atom> ligand_clusters, 
		string cluster_name, const string &names_dir, const bool for_gclus) {

		map<string, unique_ptr<Json::Value>> ligands;
		for (auto &kv : ligand_clusters) {
			int cluster_id = (kv.first == -1 ? 0 : kv.first); // unclustered is cluster #0
			Molib::Atom &atom = *kv.second;
			const Molib::Molecules &mols = atom.br().br().br().br().br().br(); // mols.name() is the nr-pdb name
			const Molib::Molecule &molecule = atom.br().br().br().br().br(); // molecule.name() is the source pdb name
			const Molib::Assembly &assembly = atom.br().br().br().br();
			const Molib::Model &model = atom.br().br().br();
			const Molib::Chain &chain = atom.br().br();
			const Molib::Residue &residue = atom.br();
			dbgmsg(mols.name());
			if (ligands.find(mols.name()) == ligands.end()) {
				ligands.insert(make_pair(mols.name(), unique_ptr<Json::Value>(new Json::Value(Json::arrayValue))));
			}
			vector<string> lig_name;
			try {
				inout::Inout::read_file(names_dir + "/" + 
				(
					(residue.rest() == Molib::Residue::protein || residue.rest() == Molib::Residue::nucleic) ? 
						(molecule.name().substr(0,4) + molecule.name().substr(5,1)) 
					: 
						residue.resn()
				), lig_name);
			}
			catch (exception& e) {
				//~ cerr << e.what() << " ... skipping ... " << endl;
				dbgmsg(e.what() << " ... skipping ... ");
			}
			ligands[mols.name()]->append(to_string(cluster_id) + ":" + residue.resn() + ":" + to_string(residue.resi()) 
				+ ":" + chain.chain_id() + ":" + to_string(model.number()) 
				+ ":" + to_string(assembly.number()) + ":" + molecule.name() + ":" + mols.name() + ":" 
				+ (lig_name.empty() ? "-" : remove_special_characters(lig_name[0])) + ":"
				+ (help::non_specific_binders.find(residue.resn()) == help::non_specific_binders.end() ? "S" : "N")
										);
			dbgmsg(residue.resn() << ":" << to_string(cluster_id) << " molecule_name = " << molecule.name() << " mols_name = " << mols.name());
		}
		for (auto &i : ligands) {
			//~ const string &pdb_id = i.first;
			const string &mols_name = i.first;
			//~ dbgmsg("pdb_id = " << pdb_id);
			dbgmsg("mols_name = " << mols_name);
			Json::Value &vec = *(i.second);
			//~ JsonReader::iterator it = jr.find({make_pair("pdb_id", pdb_id)});
			JsonReader::iterator it = jr.end();
			if(for_gclus) {
				auto vs = help::ssplit(mols_name, " ");
				 it = jr.find({make_pair("pdb_id", vs[0]), make_pair("chain_id", vs[1])});
				 dbgmsg("vs[0] = " << vs[0] << " vs[1] = " << vs[1]);
			} else {
				 it = jr.find({make_pair("pdb_id", mols_name)});
			}
			if (it != jr.end()) {
				dbgmsg("found json value " << mols_name);
				(*it)[cluster_name] = vec;
				dbgmsg("after found json value " << mols_name);
			}
		}
	}
	void remove_ligands_not_in_aligned_region(Molib::Molecules &mols, Json::Value aligned_residues) {
		set<Molib::Residue::res_tuple2> a = common_ligands::json_to_set(aligned_residues).second;
		for (auto &molecule : mols) {
			for (auto &assembly : molecule) {
				for (auto &model : assembly) {
					for (auto &chain : model) {
						// remove hetero, ion and water residues that are not close to the query
						chain.remove_if([&a, &chain, &model](const Molib::Residue &r) {
							map<Molib::Residue::res_tuple2, Molib::Residue::res_tuple2> respair;
							model.remarks(7, Molib::Residue::res_tuple2(chain.chain_id(), r.resn(), r.resi(), r.ins_code()), respair);
							set<Molib::Residue::res_tuple2> b;
							for (auto &i : respair) { b.insert(i.first); }
							dbgmsg("LIGAND " << chain.chain_id() << ":" << r.resn() << ":" << r.resi() << ":" << r.ins_code());
#ifndef NDEBUG
							for (auto &i : b) {
								dbgmsg("B " << get<0>(i) << ":" << get<1>(i) << ":" << get<2>(i) << ":" << get<3>(i));
							}
#endif
							vector<Molib::Residue::res_tuple2> v(b.size());
							vector<Molib::Residue::res_tuple2>::iterator it = set_intersection(a.begin(), a.end(), b.begin(), b.end(), v.begin());
							v.resize(it - v.begin());
#ifndef NDEBUG
							for (auto &iv : v) {
								dbgmsg("INTERSECTION " << get<0>(iv) << ":" << get<1>(iv) << ":" << get<2>(iv) << ":" << get<3>(iv));
							}
#endif
							if ((r.rest() == Molib::Residue::protein && v.size() < 5 || (r.rest() == Molib::Residue::nucleic && v.size() < 4))
								|| (r.rest() == Molib::Residue::hetero && v.size() < 3) || (r.rest() == Molib::Residue::ion && v.size() < 3)) {
								return true; // delete if less than 2 residues are in intersection
							}
							return false; // don't delete
						});
					}
				}
			}
		}
	}
	cluster::PairwiseDistances<Molib::Atom> create_pairwise_distances(
		const Molib::AtomVec &ligands, const double eps) {
		cluster::PairwiseDistances<Molib::Atom> pairwise_distances;
		Molib::MolGrid grid(ligands);
		for (auto &patom : ligands) {
			for (auto &pneighbor : grid.get_neighbors(*patom, eps)) {
				pairwise_distances[patom][pneighbor] = 
					patom->crd().distance(pneighbor->crd());
			}
		}
		return pairwise_distances;
	}
	void generate_clusters_of_ligands(const string &json_file, const string &json_with_ligs_file, 
		const string &geo_dir, const string &names_dir, const bool neighb, const double hetero_clus_rad, 
		const int hetero_min_pts, const double min_z_score, const bool for_gclus) {
		JsonReader jr;
		Molib::NRset nrset;
		//~ Molib::PDBreader pr(Molib::PDBreader::all_models);
		cluster::MapD<Molib::Atom> scores;
		// the aligned files from json
		jr.parse_JSON(json_file);
		for (auto &d : jr.root()) {
			/* if something goes wrong, e.g., a pdb file is not found, 
			 * don't exit, just skip that file..
			 */
			try { 
				const string pdb_id = d["pdb_id"].asString();
				const string chain_ids = d["chain_id"].asString();
				const double z_score = d["alignment"][0]["scores"]["z_score"].asDouble();
				if (z_score < min_z_score) continue;
				const string pdb_file = geo_dir + "/" 
					+ (for_gclus ? (pdb_id + chain_ids) : pdb_id) + ".pdb";
				dbgmsg(pdb_id + " " + chain_ids);
				//~ pr.rewind();
				Molib::PDBreader pr(pdb_file, Molib::PDBreader::all_models);
				Molib::Molecules &mols = nrset.add(new Molib::Molecules(pr.parse_molecule()));
				for (auto &molecule : mols)
					for (auto &patom : molecule.get_atoms())
						scores[patom] = z_score;
				dbgmsg("from genclus : " << mols.name());
				if (neighb) remove_ligands_not_in_aligned_region(mols, 
					d["alignment"][0]["aligned_residues"]); // remove ligands that are not near aligned residues
				mols.rotate(Geom3D::Matrix(d["alignment"][0]["rotation_matrix"], 
					d["alignment"][0]["translation_vector"]), true); // inverse rotation
			}
			catch (exception& e) {
				nrset.erase(nrset.size() - 1); // delete the last mols which is empty
				cerr << e.what() << " ... skipping ... " << endl;
			}
		}
#ifndef NDEBUG
		for (auto &mols : nrset) {
			dbgmsg("molecules name = " << mols.name());
		}
#endif
		//~ Molib::AtomVec protein_ligands = get_atoms(nrset, Molib::Residue::protein);
		//~ Molib::AtomVec nucleic_ligands = get_atoms(nrset, Molib::Residue::nucleic);
		//~ Molib::AtomVec hetero_ligands = get_atoms(nrset, Molib::Residue::hetero);
		//~ Molib::AtomVec ion_ligands = get_atoms(nrset, Molib::Residue::ion);
		Molib::AtomVec protein_ligands = nrset.get_atoms("", Molib::Residue::protein);
		Molib::AtomVec nucleic_ligands = nrset.get_atoms("", Molib::Residue::nucleic);
		Molib::AtomVec hetero_ligands = nrset.get_atoms("", Molib::Residue::hetero);
		Molib::AtomVec ion_ligands = nrset.get_atoms("", Molib::Residue::ion);
		//~ Molib::AtomVec water_ligands = get_atoms(nrset, Molib::Residue::water);
		//~ cluster::MapD<Molib::Atom> scores; // we don't use scores
		//~ cluster::MapD<Molib::Atom> scores; // we don't use scores
		cluster::PairwiseDistances<Molib::Atom> pd_protein = create_pairwise_distances(protein_ligands, 4.1);
		cluster::Optics<Molib::Atom> op(pd_protein, scores, 4.1, 3);
		cluster::PairwiseDistances<Molib::Atom> pd_nucleic = create_pairwise_distances(nucleic_ligands, 4.1);
		cluster::Optics<Molib::Atom> on(pd_nucleic, scores, 4.1, 3);
		//~ cluster::MapD<Molib::Atom> hetero_scores; for (auto &pa : hetero_ligands) hetero_scores[pa] = stod(pa->atom_name());
		cluster::PairwiseDistances<Molib::Atom> pd_hetero = create_pairwise_distances(hetero_ligands, 3.1);
		//~ cluster::Optics<Molib::Atom> oh(pd_hetero, scores, 3.1, 50);
		//~ cluster::Optics<Molib::Atom> oh(pd_hetero, scores, hetero_clus_rad + 0.1, (int) ceil(hetero_ligands.size() / hetero_min_pts_factor) + 1);
		cluster::Optics<Molib::Atom> oh(pd_hetero, scores, hetero_clus_rad + 0.1, hetero_min_pts);
		cluster::PairwiseDistances<Molib::Atom> pd_ion = create_pairwise_distances(ion_ligands, 3.1);
		cluster::Optics<Molib::Atom> oi(pd_ion, scores, 3.1, 4);
		//~ cluster::PairwiseDistances<Molib::Atom> pd_water = create_pairwise_distances(water_ligands, 3.0)
		//~ cluster::Optics<Molib::Atom> op(pd_water, scores, 3.0, 4);

		auto protein_clusters = op.extract_dbscan(4.0);
		auto nucleic_clusters = on.extract_dbscan(4.0);
		//~ auto hetero_clusters = oh.extract_dbscan(0.5);
		//~ auto hetero_clusters = oh.extract_dbscan(2.0, 100, true);
		auto hetero_clusters = oh.extract_dbscan(hetero_clus_rad, 100, true);
		auto ion_clusters = oi.extract_dbscan(3.0);
		//~ vector<pair<Molib::Atom*, int>> water_clusters = ow.extract_dbscan(2.0);
		add_to_json(jr, protein_clusters.first, "protein", names_dir, for_gclus);
		add_to_json(jr, nucleic_clusters.first, "nucleic", names_dir, for_gclus);
		add_to_json(jr, hetero_clusters.first, "hetero", names_dir, for_gclus);
		add_to_json(jr, ion_clusters.first, "ion", names_dir, for_gclus);
		//~ add_to_json(jr, water_clusters.first, "water");
		inout::Inout::file_open_put_stream(json_with_ligs_file, stringstream(jr.output_json()));
		//~ cout << jr.root() << endl;
	}
}

