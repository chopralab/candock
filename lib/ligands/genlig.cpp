#include <iostream>
#include <exception>
#include <typeinfo>
#include "jsonreader.hpp"
#include "nosqlreader.hpp"
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"
#include "geom3d/matrix.hpp"
#include "helper/error.hpp"
#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "common.hpp"
#include "genlig.hpp"
using namespace std;

ostream& operator<<(ostream& os, const genlig::BindingSiteClusters& bsites)	{
	for (auto &kv : bsites) {
		const int &cluster_number = kv.first;
		const Molib::Molecules &ligands = kv.second;
		os << "REMARK  99 ______________________ BEGINNING CLUSTER #" << cluster_number << " ______________________" << endl;
		os << ligands;
	}
	return os;
}	

namespace genlig {
	ResidueSet find_neighbor_residues(Molib::Molecule &ligand, Molib::MolGrid &grid) {
		ResidueSet neighbor_residues;
		for (auto &assembly : ligand) {
			for (auto &model : assembly) {
				for (auto &chain : model) {
					for (auto &residue : chain) {
						for (auto &atom : residue) {
							vector<Molib::Atom*> neighbors = grid.get_neighbors(atom, 4.0);
							for (auto &pa : neighbors) {
								Molib::Residue &query_residue = const_cast<Molib::Residue&>(pa->br());
								neighbor_residues.insert(&query_residue);
							}
						}
					}
				}
			}
		}
		return neighbor_residues;
	}
	ResSet find_binding_site(const Molib::Molecules &mols, const string lig_code, 
		const ResMap &aligned_to_query) {
		ResSet binding_site;
		vector<string> lc = help::ssplit(lig_code, ":");
		for (auto &molecule : mols) {
			if (molecule.name() == lc[5]) {
				for (auto &assembly : molecule) {
					if (assembly.number() == stoi(lc[4])) {
						for (auto &model : assembly) {
							if (model.number() == stoi(lc[3])) {
								ResMap respair;
								// correct: ins_code is missing from lig_code !!!
								if (model.remarks(7, Molib::Residue::res_tuple2(lc[2].at(0), 
									lc[0], stoi(lc[1]), ' '), respair)) {
									for (auto &rp : respair) {
										dbgmsg(aligned_to_query.count(rp.first));
										if (aligned_to_query.count(rp.first)) 
											binding_site.insert(aligned_to_query.at(rp.first));
									}
								}
							}
						}
					}
				}
			}
		}
		return binding_site;
	}
	Molib::Molecule get_ligand(const Molib::Molecules &mols, const string &lig_code) { // output rotated binding site residues
		Molib::Molecule ligand(lig_code);
		Molib::Assembly &lassembly = ligand.add(new Molib::Assembly(0));
		Molib::Model &lmodel = lassembly.add(new Molib::Model(1));
		char chain_id = ' ';
		Molib::Chain *plchain = nullptr;
		for (auto &molecule : mols) {
			for (auto &assembly : molecule) {
				for (auto &model : assembly) {
					for (auto &chain : model) {
						for (auto &residue : chain) {
							string lresn = (residue.rest() == Molib::Residue::protein ? 
								"PPP" : (residue.rest() == Molib::Residue::nucleic ? 
								"NNN" : residue.resn()));
							int lresi = ((residue.rest() == Molib::Residue::protein 
								|| residue.rest() == Molib::Residue::nucleic) 
								? 1 : residue.resi());
							if (lig_code.find(lresn + ":" + to_string(lresi) + ":" 
								+ chain.chain_id() + ":" + to_string(model.number()) 
								+ ":" + to_string(assembly.number()) + ":" + molecule.name() 
								+ ":" + mols.name()) != string::npos) {
									if (chain.chain_id() != chain_id) {
										plchain = &lmodel.add(new Molib::Chain(chain.chain_id()));
										chain_id = chain.chain_id();
									}
									Molib::Chain &lchain = *plchain;
									lchain.add(new Molib::Residue(residue));
									dbgmsg(lresn + ":" + to_string(lresi) + ":" + chain.chain_id() 
										+ ":" + to_string(model.number()) + ":" + to_string(assembly.number()) 
										+ ":" + molecule.name() + ":" + mols.name() + " lig_code = "
										+ lig_code);
									dbgmsg(residue);
							}
						}
					}
				}
			}
		}
		return ligand;
	}
	Molib::Molecule mark(const ResidueSet &neighbor_residues, 
		const ResSet &aligned_part_of_binding_site, const string &lig_code) {
		Molib::Molecule bsite(lig_code);
		Molib::Assembly &assembly = bsite.add(new Molib::Assembly(0));
		Molib::Model &model = assembly.add(new Molib::Model(1));
		char chain_id = ' ';
		Molib::Chain *pchain;
		for (auto &presidue : neighbor_residues) {
			if (presidue->br().chain_id() != chain_id) {
				pchain = &model.add(new Molib::Chain(presidue->br().chain_id()));
				chain_id = presidue->br().chain_id();
			}
			Molib::Chain &chain = *pchain;
			Molib::Residue &residue = chain.add(new Molib::Residue(*presidue));
			Molib::Residue::res_tuple2 r(residue.br().chain_id(), 
				help::to_string(help::get_one_letter(residue.resn())), 
				residue.resi(), residue.ins_code());
			if (aligned_part_of_binding_site.count(r))
				for (auto &atom : residue) atom.set_idatm_type("Npl");
			else
				for (auto &atom : residue) atom.set_idatm_type("C3");
			
		}
		return bsite;
	}
	Molib::MolGrid parse_receptor(const string &receptor_file, const string &receptor_chain_id) {
		// read query PDB protein
		Molib::PDBreader qpdb(receptor_file, 
			Molib::PDBreader::first_model|Molib::PDBreader::hydrogens);
		Molib::Molecules query_mols = qpdb.parse_molecule();
		return Molib::MolGrid(query_mols[0].get_atoms(receptor_chain_id, 
			Molib::Residue::protein));
	}
	void generate_ligs(Molib::MolGrid &gridrec, const string &json_file, 
		const string &bio_dir, const string &lig_code, const string &lig_file, 
		const string &bsite_file) {
		JsonReader jr;
		NosqlReader nr;
		Molib::NRset nrset;
		// read bio PDBs
		//~ Molib::PDBreader biopdb(Molib::PDBreader::all_models|Molib::PDBreader::hydrogens);
		// the aligned files from json
		jr.parse_JSON(json_file);
		for (auto &d : jr.root()) {
			try { // if something goes wrong, e.g., pdb file is not found, don't exit..
				const string pdb_id = d["pdb_id"].asString();
				const string chain_ids = d["chain_id"].asString();
				//~ const string mols_name = pdb_id + " " + chain_ids;
				//~ const string pdb_file = bio_dir + "/" + pdb_id + chain_ids + ".pdb";
				const string mols_name = pdb_id;
				const string pdb_file = bio_dir + "/" + pdb_id + ".pdb";
				if (lig_code.find(mols_name) != string::npos) {
					//~ Molib::Molecules &mols = nrset.add(new Molib::Molecules());
					//~ biopdb.parse_PDB(mols, pdb_file);
					//~ biopdb.rewind();
					Molib::PDBreader biopdb(pdb_file, 
						Molib::PDBreader::all_models|Molib::PDBreader::hydrogens);
					Molib::Molecules &mols = 
						nrset.add(new Molib::Molecules(biopdb.parse_molecule()));
					dbgmsg(Geom3D::Matrix(d["alignment"][0]["rotation_matrix"], 
						d["alignment"][0]["translation_vector"]));
					dbgmsg(mols_name);
					mols.rotate(Geom3D::Matrix(d["alignment"][0]["rotation_matrix"], 
						d["alignment"][0]["translation_vector"]), true); // inverse rotation
					ResMap aligned_to_query = common_ligands::json_to_map_reverse(d["alignment"][0]["aligned_residues"]);
#ifndef NDEBUG
					for (auto &kv : aligned_to_query) {
						auto &aligned = kv.first;							
						auto &query = kv.second;
						dbgmsg("set " << get<0>(query) << ":" << get<1>(query) 
							<< ":" << get<2>(query) << ":" << get<3>(query));
						dbgmsg("set " << get<0>(aligned) << ":" << get<1>(aligned) 
							<< ":" << get<2>(aligned) << ":" << get<3>(aligned));
					}
#endif
					ResSet aligned_part_of_binding_site = find_binding_site(mols, lig_code, aligned_to_query);
					Molib::Molecules ligands;
					ligands.add(new Molib::Molecule(get_ligand(mols, lig_code)));
					ResidueSet neighbor_residues = find_neighbor_residues(ligands.first(), gridrec);
					Molib::Molecules binding_sites;
					binding_sites.add(new Molib::Molecule(mark(neighbor_residues, 
						aligned_part_of_binding_site, lig_code)));
					inout::output_file(ligands, lig_file);
					inout::output_file(binding_sites, bsite_file);
					break; // the end
				}
			}
			catch (exception& e) {
				nrset.erase(nrset.size() - 1); // delete the last mols which is empty
				cerr << e.what() << " ... skipping ... " << endl;
			}
		}
	}
	void generate_ligands(const string &receptor_file, const string &receptor_chain_id, 
		const string &json_file, const string &bio_dir, const string &lig_code,
		const string &lig_file, const string &bsite_file) {
		Molib::MolGrid gridrec = parse_receptor(receptor_file, receptor_chain_id);
		generate_ligs(gridrec, json_file, bio_dir, lig_code, lig_file, bsite_file);
	
	}
	BindingSiteClusters generate_binding_site_prediction(const string &json_with_ligs_file, 
		const string &bio_dir, const int num_bsites) {
		JsonReader jr;
		BindingSiteClusters bsites;
		dbgmsg("starting generate binding site prediction");
		// read hetero ligands from json_with_ligs file
		jr.parse_JSON(json_with_ligs_file);
		map<string, Geom3D::Matrix> bio_file_to_matrix;
		map<string, vector<pair<int, string>>> ligands_by_bio_file;
		for (auto &d : jr.root()) {
			dbgmsg(Geom3D::Matrix(d["alignment"][0]["rotation_matrix"], 
				d["alignment"][0]["translation_vector"]));
			Geom3D::Matrix mx(d["alignment"][0]["rotation_matrix"], 
				d["alignment"][0]["translation_vector"]);
			const Json::Value &hetero = d["hetero"];
			for (int i = 0; i < hetero.size(); ++i) { // go over the hetero ligands only
				dbgmsg(hetero[i].asString());
				const vector<string> ligand = help::ssplit(hetero[i].asString(), ":");
				dbgmsg(ligand.size());
				const int cluster_number = stoi(ligand[0]);
				const string resn = ligand[1];
				const int resi = stoi(ligand[2]);
				const string chain_id = ligand[3];
				const int model_number = stoi(ligand[4]);
				const int assembly_number = stoi(ligand[5]);
				const string molecule_name = ligand[6];
				const string mols_name = ligand[7];
				const string name = ligand[8];

				const string bio_file = bio_dir + "/" + mols_name + ".pdb";
				// lig_code = resn:resi:chain_id:model_num:assembly_num:molecule_name
				const string lig_code = resn + ":" + help::to_string(resi)
					+ ":" + chain_id + ":" + help::to_string(model_number)
					+ ":" + help::to_string(assembly_number)
					+ ":" + molecule_name +  ":" + mols_name;
				dbgmsg("after lig_code");
				
				ligands_by_bio_file[bio_file].push_back({cluster_number, lig_code});
				bio_file_to_matrix[bio_file] = mx;
			}
		}
		// read bio PDBs
		//~ Molib::PDBreader biopdb(Molib::PDBreader::all_models);
		for (auto &kv1 : ligands_by_bio_file) {
			const string &bio_file = kv1.first;
			const Geom3D::Matrix &mx = bio_file_to_matrix[bio_file];
			try { // if something goes wrong, e.g., pdb file is not found, don't exit..
				//~ biopdb.rewind();
				Molib::PDBreader biopdb(bio_file, 
					Molib::PDBreader::all_models);
				Molib::Molecules mols = biopdb.parse_molecule();
				mols.rotate(mx, true); // inverse rotation
				for (auto &kv2 : kv1.second) {
					const int cluster_number = kv2.first;
					if (cluster_number <= num_bsites) { // trim number of binding sites
						const string &lig_code = kv2.second;
						bsites[cluster_number].add(new Molib::Molecule(get_ligand(mols, lig_code)));
					}
				}
			}
			catch (exception& e) {
				cerr << e.what() << " ... skipping ... " << endl;
			}
		}
		return bsites;
	}
}
