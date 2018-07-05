#include <iostream>
#include <exception>
#include <typeinfo>
#include "candock/molib/grid.hpp"
#include "candock/ligands/jsonreader.hpp"
#include "candock/ligands/nosqlreader.hpp"
#include "candock/molib/nrset.hpp"
#include "candock/parser/fileparser.hpp"
#include "candock/helper/help.hpp"
#include "candock/geometry/matrix.hpp"
#include "candock/helper/error.hpp"
#include "candock/helper/inout.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/path.hpp"
#include "candock/ligands/common.hpp"
#include "candock/ligands/genlig.hpp"
using namespace std;

namespace candock{
namespace genlig {
	ResidueSet find_neighbor_residues(molib::Molecule &ligand, molib::Atom::Grid &grid) {
		ResidueSet neighbor_residues;
		for (auto &assembly : ligand) {
			for (auto &model : assembly) {
				for (auto &chain : model) {
					for (auto &residue : chain) {
						for (auto &atom : residue) {
							dbgmsg("before get neighbors");
							molib::Atom::Vec neighbors = grid.get_neighbors(atom.crd(), 4.0);
							dbgmsg("number of neighbors = " << neighbors.size());
							for (auto &pa : neighbors) {
								dbgmsg("pa = " << *pa);
								molib::Residue &query_residue = const_cast<molib::Residue&>(pa->br());
								dbgmsg("query_residue = " << query_residue);
								neighbor_residues.insert(&query_residue);
							}
						}
					}
				}
			}
		}
		return neighbor_residues;
	}
	ResSet find_binding_site(const molib::Molecules &mols, const string lig_code, 
		const ResMap &aligned_to_query) {
		ResSet binding_site;
		vector<string> lc = help::ssplit(lig_code, ":");
		for (auto &molecule : mols) {
			dbgmsg("molecule.name = " << molecule.name() << " lc[5[ = " << lc[5]);
			if (molecule.name() == lc[5]) {
				for (auto &assembly : molecule) {
					if (assembly.number() == stoi(lc[4])) {
						for (auto &model : assembly) {
							if (model.number() == stoi(lc[3])) {
								ResMap respair;
								// correct: ins_code is missing from lig_code !!!
								if (model.remarks(7, molib::Residue::res_tuple2(lc[2].at(0), 
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
	molib::Molecule get_ligand(const molib::Molecules &mols, const string &lig_code) { // output rotated binding site residues
		molib::Molecule ligand(lig_code);
		molib::Assembly &lassembly = ligand.add(new molib::Assembly(0));
		molib::Model &lmodel = lassembly.add(new molib::Model(1));
		char chain_id = ' ';
		molib::Chain *plchain = nullptr;
		for (auto &molecule : mols) {
			for (auto &assembly : molecule) {
				for (auto &model : assembly) {
					for (auto &chain : model) {
						for (auto &residue : chain) {
							string lresn = (residue.rest() == molib::Residue::protein ? 
								"PPP" : (residue.rest() == molib::Residue::nucleic ? 
								"NNN" : residue.resn()));
							int lresi = ((residue.rest() == molib::Residue::protein 
								|| residue.rest() == molib::Residue::nucleic) 
								? 1 : residue.resi());
							if (lig_code.find(lresn + ":" + to_string(lresi) + ":" 
								+ chain.chain_id() + ":" + to_string(model.number()) 
								+ ":" + to_string(assembly.number()) + ":" + molecule.name() 
								+ ":" + mols.name()) != string::npos) {
									if (chain.chain_id() != chain_id) {
										plchain = &lmodel.add(new molib::Chain(chain.chain_id()));
										chain_id = chain.chain_id();
									}
									molib::Chain &lchain = *plchain;
									lchain.add(new molib::Residue(residue));
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
	molib::Molecule mark(const ResidueSet &neighbor_residues, 
		const ResSet &aligned_part_of_binding_site, const string &lig_code) {
		molib::Molecule bsite(lig_code);
		molib::Assembly &assembly = bsite.add(new molib::Assembly(0));
		molib::Model &model = assembly.add(new molib::Model(1));
		char chain_id = ' ';
		molib::Chain *pchain = nullptr;
		for (auto &presidue : neighbor_residues) {
			if (presidue->br().chain_id() != chain_id) {
				pchain = &model.add(new molib::Chain(presidue->br().chain_id()));
				chain_id = presidue->br().chain_id();
			}
			molib::Chain &chain = *pchain;
			molib::Residue &residue = chain.add(new molib::Residue(*presidue));
			molib::Residue::res_tuple2 r(residue.br().chain_id(), 
				std::to_string(help::get_one_letter(residue.resn())), 
				residue.resi(), residue.ins_code());
			if (aligned_part_of_binding_site.count(r))
				for (auto &atom : residue) atom.set_idatm_type("Npl");
			else
				for (auto &atom : residue) atom.set_idatm_type("C3");
			
		}
		return bsite;
	}
	//~ molib::Atom::Grid parse_receptor(const string &receptor_file, const string &receptor_chain_id) {
		//~ // read query PDB protein
		//~ molib::PDBreader qpdb(receptor_file, 
			//~ molib::PDBreader::first_model|molib::PDBreader::hydrogens);
		//~ molib::Molecules query_mols = qpdb.parse_molecule();
		//~ return molib::Atom::Grid(query_mols[0].get_atoms(receptor_chain_id, 
			//~ molib::Residue::protein));
	//~ }
	void generate_ligands(const string &receptor_file, const string &receptor_chain_id, 
		const string &json_file, const string &bio_dir, const string &lig_code,
		const string &lig_file, const string &bsite_file) {
	//~ void generate_ligs(molib::Atom::Grid &gridrec, const string &json_file, 
		//~ const string &bio_dir, const string &lig_code, const string &lig_file, 
		//~ const string &bsite_file) {

		// read query PDB protein
		parser::FileParser qpdb(receptor_file, 
			parser::first_model|parser::hydrogens);
		molib::Molecules query_mols = qpdb.parse_molecule();
		molib::Atom::Grid gridrec(query_mols[0].get_atoms(receptor_chain_id, 
			molib::Residue::protein));

		JsonReader jr;
		molib::NRset nrset;
		// read bio PDBs
		//~ molib::PDBreader biopdb(molib::PDBreader::all_models|molib::PDBreader::hydrogens);
		// the aligned files from json
		jr.parse_JSON(json_file);
		for (auto &d : jr.root()) {
			try { // if something goes wrong, e.g., pdb file is not found, don't exit..
				const string pdb_id = d["pdb_id"].asString();
				const string chain_ids = d["chain_id"].asString();
				const string mols_name = pdb_id + chain_ids;
				const string pdb_file = Path::join(bio_dir, pdb_id + chain_ids + ".pdb");
				//~ const string mols_name = pdb_id;
				//~ const string pdb_file = bio_dir + "/" + pdb_id + ".pdb";
				//~ cout << mols_name << " " << pdb_file << endl;
				if (lig_code.find(mols_name) != string::npos) {
					//~ molib::Molecules &mols = nrset.add(new molib::Molecules());
					//~ biopdb.parse_PDB(mols, pdb_file);
					//~ biopdb.rewind();
					parser::FileParser biopdb(pdb_file, 
						parser::all_models|parser::hydrogens);
					molib::Molecules &mols = 
						nrset.add(new molib::Molecules(biopdb.parse_molecule()));
					dbgmsg(geometry::Matrix(d["alignment"][0]["rotation_matrix"], 
						d["alignment"][0]["translation_vector"]));
					dbgmsg(mols_name);
					mols.rotate(geometry::Matrix(d["alignment"][0]["rotation_matrix"], 
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
					molib::Molecules ligands;
					ligands.add(new molib::Molecule(get_ligand(mols, lig_code)));
					dbgmsg(ligands);
					ResidueSet neighbor_residues = find_neighbor_residues(ligands.first(), gridrec);
					dbgmsg(neighbor_residues);
					molib::Molecules binding_sites;
					binding_sites.add(new molib::Molecule(mark(neighbor_residues, 
						aligned_part_of_binding_site, lig_code)));
					Inout::output_file(ligands, lig_file);
					Inout::output_file(binding_sites, bsite_file);
					break; // the end
				}
			}
			catch (exception& e) {
				nrset.erase(nrset.size() - 1); // delete the last mols which is empty
				log_warning << e.what() << " ... skipping ... " << endl;
			}
		}
	}
	//~ void generate_ligands(const string &receptor_file, const string &receptor_chain_id, 
		//~ const string &json_file, const string &bio_dir, const string &lig_code,
		//~ const string &lig_file, const string &bsite_file) {
//~ 
//~ 
		//~ generate_ligs(gridrec, json_file, bio_dir, lig_code, lig_file, bsite_file);
	//~ 
	//~ }
	pair<BindingSiteClusters, BindingSiteScores> generate_binding_site_prediction(const string &json_with_ligs_file, 
		const string &bio_dir, const int num_bsites) {
		JsonReader jr;
		BindingSiteClusters bsites;
		BindingSiteScores bscores;
		dbgmsg("starting generate binding site prediction");
		// read hetero ligands from json_with_ligs file
		jr.parse_JSON(json_with_ligs_file);
		map<string, geometry::Matrix> bio_file_to_matrix;
		struct LigandInfo {
			int cluster_number;
			string lig_code;
			double z_score;
		};
		map<string, vector<LigandInfo>> ligands_by_bio_file;
		for (auto &d : jr.root()) {
			dbgmsg(geometry::Matrix(d["alignment"][0]["rotation_matrix"], 
				d["alignment"][0]["translation_vector"]));
			geometry::Matrix mx(d["alignment"][0]["rotation_matrix"], 
				d["alignment"][0]["translation_vector"]);
			const double z_score = d["alignment"][0]["scores"]["z_score"].asDouble();
			const Json::Value &hetero = d["hetero"];
			for (size_t i = 0; i < hetero.size(); ++i) { // go over the hetero ligands only
				dbgmsg(hetero[static_cast<int>(i)].asString());
				const vector<string> ligand = help::ssplit(hetero[ static_cast<int>(i)].asString(), ":");
			//	dbgmsg(ligand.size());
				const int cluster_number = stoi(ligand[0]);
				const string resn = ligand[1];
				const int resi = stoi(ligand[2]);
				const string chain_id = ligand[3];
				const int model_number = stoi(ligand[4]);
				const int assembly_number = stoi(ligand[5]);
				const string molecule_name = ligand[6];
				const string mols_name = ligand[7];
				const string name = ligand[8];

				const string bio_file = Path::join(bio_dir, mols_name + ".pdb");
				// lig_code = resn:resi:chain_id:model_num:assembly_num:molecule_name
				const string lig_code = resn + ":" + std::to_string(resi)
					+ ":" + chain_id + ":" + std::to_string(model_number)
					+ ":" + std::to_string(assembly_number)
					+ ":" + molecule_name +  ":" + mols_name;
				dbgmsg("after lig_code");
				
				ligands_by_bio_file[bio_file].push_back(LigandInfo{cluster_number, lig_code, z_score});
				bio_file_to_matrix[bio_file] = mx;
			}
		}
		// read bio PDBs
		//~ molib::PDBreader biopdb(molib::PDBreader::all_models);
		for (auto &kv1 : ligands_by_bio_file) {
			const string &bio_file = kv1.first;
			const geometry::Matrix &mx = bio_file_to_matrix[bio_file];
			try { // if something goes wrong, e.g., pdb file is not found, don't exit..
				//~ biopdb.rewind();
				parser::FileParser biopdb(bio_file, 
					parser::all_models);
				molib::Molecules mols = biopdb.parse_molecule();
				mols.rotate(mx, true); // inverse rotation
				for (auto &linf : kv1.second) {
					if (linf.cluster_number <= num_bsites) { // trim number of binding sites
						bsites[linf.cluster_number].add(new molib::Molecule(get_ligand(mols, linf.lig_code)));
						bscores[linf.cluster_number] = max(bscores[linf.cluster_number], linf.z_score);
					}
				}
			}
			catch (exception& e) {
				log_warning << e.what() << " ... skipping ... " << endl;
			}
		}
		return {bsites, bscores};
	}
}

namespace molib {
	ostream& operator<<(ostream& os, const map<int, molib::Molecules>& bsites) {
		for (auto &kv : bsites) {
			const int &cluster_number = kv.first;
			const molib::Molecules &ligands = kv.second;
			os << "REMARK  99 ______________________ BEGINNING CLUSTER #" << cluster_number << " ______________________" << endl;
			os << ligands;
		}
		return os;
	}
}

ostream& operator<<(ostream& os, const map<int, double>& bscores) {
	for (auto &kv : bscores) {
		const int &cluster_number = kv.first;
		const double &z_score = kv.second;
		os << cluster_number << " " << z_score << endl;
	}
	return os;
}
}
