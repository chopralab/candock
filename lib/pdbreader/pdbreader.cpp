#include <memory>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "geom3d/coordinate.hpp"
#include "helper/inout.hpp"
#include "helper/help.hpp"
#include "helper/debug.hpp"
#include "molecule.hpp"
#include "pdbreader.hpp"
#include "bond.hpp"
#include "fragmenter/fragmenter.hpp"
using namespace std;

namespace Molib {
	void PDBreader::rewind() { p->set_pos(0); }
	void PDBreader::set_flags(unsigned int hm) { p->set_hm(hm); }
	void PDBreader::parse_molecule(Molecules &mols) { 
		p->parse_molecule(mols); 
		dbgmsg("PARSED MOLECULES : " << endl << mols);
	}
	Molecules PDBreader::parse_molecule() { 
		Molecules mols; 
		p->parse_molecule(mols); 
		dbgmsg("PARSED MOLECULES : " << endl << mols);
		return mols; 
	}
	
	PDBreader::PDBreader(const string &molecule_file, unsigned int hm, 
		const int num_occur) { 
		auto ret = molecule_file.find_last_of(".");
		if (ret == string::npos) 
			throw Error("die : could not determine the file type of the input molecule");
		string extension = molecule_file.substr(ret + 1);
		transform(extension.begin(), extension.end(), 
			extension.begin(), ::toupper);
		dbgmsg("pdb reader options = " << hm << " molecule is a " << boolalpha 
			<< (extension == "PDB" || extension == "ENT" ? "PDB file" :
			(extension == "MOL2" ? "Mol2 file" : "undetermined file type")));
		if (extension == "PDB" || extension == "ENT")
			p = new PdbParser(molecule_file, hm, num_occur);
		else if (extension == "MOL2")
			p = new Mol2Parser(molecule_file, hm, num_occur);
		else
			throw Error("die : could not determine the file type of the input molecule");
	}
	void PDBreader::Parser::__generate_molecule(Molecules &mols, bool &found_molecule, const string &name) {
		if(!found_molecule) { // if there were no REMARK or BIOMOLECULE...
			mols.add(new Molecule(name));
			found_molecule = true;
		}
	}
	void PDBreader::Parser::__generate_assembly(Molecules &mols, bool &found_assembly, int assembly_number, const string &name) {
		if(!found_assembly) { // if there were no REMARK or BIOMOLECULE...
			mols.last().add(new Assembly(assembly_number, name));
			found_assembly = true;
		}
	}	
	void PDBreader::Parser::__generate_model(Molecules &mols, bool &found_model, int model_number) {
		if(!found_model) { // if there were no REMARK or BIOMOLECULE...
			mols.last().last().add(new Model(model_number));
			found_model = true;
		}
	}
	void PDBreader::Mol2Parser::parse_molecule(Molecules &mols) {
		vector<string> mol2_raw;
		__pos = inout::Inout::read_file(__molecule_file, mol2_raw, 
			inout::Inout::panic, __pos, __num_occur, "@<TRIPOS>MOLECULE");
		bool found_molecule = false, found_assembly = false, found_model = false;
		map<const Model*, map<const int, Atom*>> atom_number_to_atom;
		for (int i = 0; i < mol2_raw.size(); ++i) {
			const string &line = mol2_raw[i];
			if (line.find("@<TRIPOS>ATOM") != string::npos) {
				__generate_molecule(mols, found_molecule, "");
				__generate_assembly(mols, found_assembly, 0, "ASYMMETRIC UNIT");
				__generate_model(mols, found_model, 1);
				for (i = i + 1; i < mol2_raw.size(); ++i) {
					const string &atom_line = mol2_raw[i];
					dbgmsg(atom_line);
					if (atom_line.empty() || atom_line[0] == '@') break;
					stringstream ss(atom_line);
					int atom_id, subst_id;
					string atom_name, atom_type, subst_name;
					double x, y, z, charge;
					ss >> atom_id >> atom_name >> x >> y >> z >> atom_type 
						>> subst_id >> subst_name >> charge;
					dbgmsg("atom_id = " << atom_id << " atom_name = " << atom_name
						<< " x = " << x << " y = " << y << " z = " << z
						<< " atom_type = " << atom_type << " subst_id = " << subst_id
						<< " subst_name = " << subst_name << " charge = " << charge);
					Geom3D::Coordinate crd(x, y, z);
					dbgmsg("is hydrogen = " << (__hm & PDBreader::hydrogens) 
						<< " atom_name = " << atom_name.at(0));
					if ((__hm & PDBreader::hydrogens) || atom_name.at(0) != 'H') {					
						Residue::res_type rest(Residue::hetero);
						if (!__giant_molecule || rest != Residue::protein || atom_name == "CA") {
							Model &model = mols.last().last().last();
							Chain *chain;
							if (!model.has_chain('A'))
								chain = &model.add(new Chain('A'));
							Residue *residue;
							if (!chain->has_residue(Residue::res_pair(subst_id, ' ')))
								residue = &chain->add(new Residue(subst_name, subst_id, ' ', rest));
							Atom &a = residue->add(new Atom(atom_id, 
															atom_name, 
															crd, 
															help::idatm_mask.at("???")
														));
							atom_number_to_atom[&model][atom_id] = &a;
						}
					}
				}
				--i;
			}
			else if(line.find("@<TRIPOS>MOLECULE") != string::npos) { 
				if (i + 1 >= mol2_raw.size()) throw Error("die : wrong format of mol2 file");
				const string &name = mol2_raw[i + 1];
				dbgmsg(name);
				found_molecule = false;
				found_assembly = false;
				found_model = false;
				__generate_molecule(mols, found_molecule, name);
				__generate_assembly(mols, found_assembly, 0, "ASYMMETRIC UNIT");
				__generate_model(mols, found_model, 1);
			}
			else if (line.find("@<TRIPOS>BOND") != string::npos) {
				Molecule &molecule = mols.last(); // the CONECT words are valid for the whole molecule
				for (i = i + 1; i < mol2_raw.size(); ++i) {
					const string &bond_line = mol2_raw[i];
					if (bond_line.empty() || bond_line[0] == '@') break;
					int bond_id, origin_atom_id, target_atom_id;
					string bond_type;
					stringstream ss(bond_line);
					ss >> bond_id >> origin_atom_id >> target_atom_id >> bond_type;
					for (auto &assembly : molecule) {
						for (auto &model : assembly) {
							auto it1 = atom_number_to_atom[&model].find(origin_atom_id);
							auto it2 = atom_number_to_atom[&model].find(target_atom_id);
							if (it1 != atom_number_to_atom[&model].end()
								&& it2 != atom_number_to_atom[&model].end()) {
								Atom &a1 = *it1->second;
								Atom &a2 = *it2->second;
								a1.connect(a2);
							}
						}
					}
				}
				connect_bonds(get_bonds_in(molecule.get_atoms()));
				--i;
			}
		}
	}
	void PDBreader::PdbParser::parse_molecule(Molecules &mols) {
		vector<string> pdb_raw;
		__pos = inout::Inout::read_file(__molecule_file, pdb_raw, inout::Inout::panic, 
			__pos, __num_occur, "REMARK   5 END");
		boost::smatch m;
		set<char> bio_chain;
		int biomolecule_number = -1;
		int coord_index = 1;
		bool found_molecule = false, found_assembly = false, found_model = false;
		map<const Model*, map<const int, Atom*>> atom_number_to_atom;
		for (string &line : pdb_raw) {
			if ((__hm & PDBreader::skip_hetatm) && line.compare(0, 6, "HETATM") == 0) continue;
			if (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0) {
				__generate_molecule(mols, found_molecule, "");
				__generate_assembly(mols, found_assembly, 0, "ASYMMETRIC UNIT");
				__generate_model(mols, found_model, 1);
				string record_name = boost::algorithm::trim_copy(line.substr(0,6));
				int atom_number = stoi(line.substr(6, 5));
				string atom_name = boost::algorithm::trim_copy(line.substr(12,4));
				char alt_loc = line.at(16);
				string resn = boost::algorithm::trim_copy(line.substr(17,3));
				char chain_id = line.at(21);
				int resi = stoi(line.substr( 22, 4));
				char ins_code = line.at(26);			
				double x_coord = stof(line.substr(30, 8));
				double y_coord = stof(line.substr(38, 8));
				double z_coord = stof(line.substr(46, 8));
				Geom3D::Coordinate crd(x_coord, y_coord, z_coord);
				string element = line.size() > 77 ? boost::algorithm::trim_copy(line.substr(76,2)) : "";
				string idatm_type = line.size() > 80 ? boost::algorithm::trim_copy(line.substr(80, 5)) : "???";
				idatm_type = help::idatm_mask.count(idatm_type) ? idatm_type : "???";
				string gaff_type = line.size() > 85 ? boost::algorithm::trim_copy(line.substr(85, 5)) : "???";
				string rest_of_line = line.size() > 90 ? boost::algorithm::trim_copy(line.substr(90)) : "";
				dbgmsg("is hydrogen = " << (__hm & PDBreader::hydrogens) << " atom_name = " << atom_name.at(0));
				if ((__hm & PDBreader::hydrogens) || atom_name.at(0) != 'H') {
					if (alt_loc == ' ' || alt_loc == 'A') {
						if (!mols.last().is_modified(Residue::res_tuple2(chain_id, resn, resi, ins_code))) {
							Residue::res_type rest;
							if (help::amino_acids.find(resn) != help::amino_acids.end()) rest = Residue::protein;
							else if (help::nucleic_acids.find(resn) != help::nucleic_acids.end()) rest = Residue::nucleic;
							else if (help::ions.find(resn) != help::ions.end()) rest = Residue::ion;
							else if (resn == "HOH") rest = Residue::water;
							else rest = Residue::hetero;  // if nothing known then hetero...
							if (!__giant_molecule || rest != Residue::protein || atom_name == "CA") {
								Model &model = mols.last().last().last();
								Chain *chain;
								if (!model.has_chain(chain_id))
									chain = &model.add(new Chain(chain_id));
								Residue *residue;
								if (!chain->has_residue(Residue::res_pair(resi, ins_code)))
									residue = &chain->add(new Residue(resn, resi, ins_code, rest));
								Atom &a = residue->add(new Atom(atom_number, 
																atom_name, 
																crd, 
																help::idatm_mask.at(idatm_type),
																element
															));
								if (gaff_type != "???")
									a.set_gaff_type(gaff_type);
								if (rest_of_line != "")
									a.set_members(rest_of_line);
								atom_number_to_atom[&model][atom_number] = &a;
							}
						}
					}
				}
			}
			else if(line.compare(0, 14, "REMARK   4 NRP") == 0) { 
				string name = "";
				if (boost::regex_search(line, m, boost::regex("REMARK   4 NRPDB\\s+(.*)"))) {
					if (m[1].matched) {
						name = m[1].str();
						dbgmsg(name);
					}
				}
				mols.set_name(name);
			}
			else if(line.compare(0, 14, "REMARK   5 MOL") == 0) { 
				string name = "";
				if (boost::regex_search(line, m, boost::regex("REMARK   5 MOLECULE\\s+(.*)"))) {
					if (m[1].matched) {
						name = m[1].str();
						dbgmsg(name);
					}
				}
				found_molecule = false;
				__generate_molecule(mols, found_molecule, name);
			}
			else if(line.compare(0, 14, "REMARK   6 BIO") == 0 || line.compare(0, 14, "REMARK   6 ASY") == 0) {
				int number = 0;
				string name = "";
				if (boost::regex_search(line, m, boost::regex("REMARK   6 (\\S+ \\S+)\\s+(\\d+)"))) {
					if (m[1].matched && m[2].matched) {
						name = m[1].str();
						number = stoi(m[2].str());
					}
				}
				found_assembly = false;
				__generate_assembly(mols, found_assembly, number, name);
			}
			else if (line.compare(0, 5, "MODEL") == 0) {
				if (boost::regex_search(line, m, boost::regex("MODEL\\s+(\\d+)"))) {
					if (m[1].matched) {
						__generate_molecule(mols, found_molecule, "");
						__generate_assembly(mols, found_assembly, 0, "ASYMMETRIC UNIT");
						if (__hm & PDBreader::all_models) found_model = false;
						__generate_model(mols, found_model, stoi(m[1].str()));
					}
				} else throw Error("die : model number not found in pdb");
			}
			else if (line.compare(0, 6, "ENDMDL") == 0) {
				if (__hm & PDBreader::first_model) return;
			}
			else if(line.compare(0, 10, "REMARK 350") == 0) { 
				if (boost::regex_search(line, m, boost::regex("REMARK 350 BIOMOLECULE:\\s*(\\d+)"))) {
					if (m[1].matched ) {
						biomolecule_number = stoi(m[1].str());
						if(biomolecule_number == 1) { // if there were no REMARKs...
							__generate_molecule(mols, found_molecule, "");
						}
					}
					bio_chain.clear();
				}
				else if (boost::regex_search(line, m, boost::regex("REMARK 350 APPLY THE FOLLOWING TO CHAINS:\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?"))
						|| boost::regex_search(line, m, boost::regex("REMARK 350                    AND CHAINS:\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?"))) {
					if (m.empty()) { throw Error("die : problem in REMARK 350 APPLY"); }
					for (auto &el : m) {
						if (el != *m.begin() && el.matched) {
							bio_chain.insert(el.str().at(0));
						}
					}
				}
				else if (boost::regex_search(line, m, boost::regex("REMARK 350   BIOMT(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)"))) {
					if (m.empty()) { throw Error("die : problem in REMARK 350   BIOMT"); }
					if (m[1].matched && m[2].matched && m[3].matched && m[4].matched && m[5].matched && m[6].matched) {
						if (biomolecule_number == -1) { throw Error("die : missing 'REMARK 350 BIOMOLECULE'"); }
						int row = stoi(m[1].str());
						int rn = stoi(m[2].str());
						double a = stod(m[3].str());
						double b = stod(m[4].str());
						double c = stod(m[5].str());
						double w = stod(m[6].str());
						if (rn < 10) {
							if (row == 1) {
								mols.last().add_bio_chain(biomolecule_number, bio_chain);
							}
							mols.last().add_bio_rota(biomolecule_number, rn, row - 1, make_tuple(a, b, c, w));
						} else {
							__giant_molecule = true; // if giant molecule, read only CA
						}
					}
				}
			}
			else if (line.compare(0, 6, "MODRES") == 0) {
				string resn = line.substr(24, 3);
				if (help::amino_acids.find(resn) != help::amino_acids.end()) {
					string resn_mod = line.substr(12, 3);
					char chain_id = line.at(16);
					int resi = stoi(line.substr(18, 4));
					char ins_code = line.at(22);
					__generate_molecule(mols, found_molecule, "");
					mols.last().add_modified(Residue::res_tuple2(chain_id, resn_mod, resi, ins_code)); // modified 'standard' residue gets saved,e.g., MET --> MSE, glycosylated...
				}
			}
			else if (line.compare(0, 4, "SITE") == 0 && boost::regex_search(line, m, boost::regex("SITE\\s+\\d+\\s+(\\S+)\\s+\\S+(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?"))) {
				if (m.empty()) { throw Error("die : problem in SITE"); }
				string site_name = (m[1].matched ? m[1].str() : "");
				for (auto &el : m) {
					if (el != *m.begin() && el.matched) {
						string s = boost::algorithm::trim_copy(el.str());
						if (s.size() > 8) {
							string resn = s.substr(0, 3);
							char chain_id = s.at(4);
							int resi = stoi(s.substr(5, 4));
							char ins_code =  (s.size() == 10 ? s.at(9) : ' ');
							if (help::amino_acids.find(resn) != help::amino_acids.end()) {
								__generate_molecule(mols, found_molecule, "");
								mols.last().add_site(site_name, Residue::res_tuple2(chain_id, resn, resi, ins_code));
							}
						}
					}
				}
			}
			else if (line.compare(0, 10, "REMARK   7") == 0) {
				string text(line.substr(11));
				vector<string> vs = help::ssplit(text, ",");
				if (vs.empty()) throw("die: remark 7 is empty");
				vector<string> token = help::ssplit(vs[0], ":");
				char lchain_id = token[0].at(0);
				string lresn = token[1];
				int lresi = stoi(token[2]);  
				char lins_code = token[3].at(0);
				Model &model = mols.last().last().last();
				vector<string>::iterator it = vs.begin();
				it++;
				for ( ; it != vs.end(); it++) {
					string &s = *it;
					vector<string> tok = help::ssplit(s, "=");
					vector<string> tok1 = help::ssplit(tok[0], ":");
					vector<string> tok2 = help::ssplit(tok[1], ":");
					dbgmsg("ADDING TO REMARK " << "[" << lchain_id << ":" << lresn << ":" << lresi << ":" << lins_code << "] => [" 
						<< tok1[0].at(0) << ":" << tok1[1] << ":" << stoi(tok1[2]) << ":" << tok1[3].at(0)
						<< "="
						<< tok2[0].at(0) << ":" << tok2[1] << ":" << stoi(tok2[2]) << ":" << tok2[3].at(0)
						<< "]");
					model.set_remark(7, Residue::res_tuple2(lchain_id, lresn, lresi, lins_code), 
							 make_pair(Residue::res_tuple2(tok1[0].at(0), tok1[1], stoi(tok1[2]), tok1[3].at(0)), 
										Residue::res_tuple2(tok2[0].at(0), tok2[1], stoi(tok2[2]), tok2[3].at(0))));
				}
			}
			else if (line.compare(0, 16, "REMARK   8 RIGID") == 0) {
				Model &model = mols.last().last().last();
				vector<string> vs = help::ssplit(line.substr(17), " ");
				AtomSet as;
				for (int i = 1; i < vs.size(); i++)	as.insert(atom_number_to_atom[&model][stoi(vs[i])]);
				model.get_rigid()[vs[0]].insert(as);
				dbgmsg(vs[0]);
#ifndef NDEBUG
				for (int i = 1; i < vs.size(); i++) dbgmsg(vs[i]);
				for (auto &a : as) dbgmsg(*a);
				dbgmsg(model.get_rigid());
#endif
			}
			else if (line.compare(0, 15, "REMARK   8 SEED") == 0) {
				Model &model = mols.last().last().last();
				vector<string> vs = help::ssplit(line.substr(16), " ");
				AtomSet as;
				for (int i = 1; i < vs.size(); i++)	as.insert(atom_number_to_atom[&model][stoi(vs[i])]);
				model.get_seeds()[vs[0]].insert(as);
				dbgmsg(model.get_seeds());
			}
			else if (line.compare(0, 15, "REMARK   8 ROTA") == 0) {
				Model &model = mols.last().last().last();
				vector<string> vs = help::ssplit(line.substr(16), " ");
				Atom &a1 = *atom_number_to_atom[&model][stoi(vs[1])];
				Atom &a2 = *atom_number_to_atom[&model][stoi(vs[2])];
				a1.connect(a2).set_rotatable(vs[0]);
			}
			else if (line.compare(0, 19, "REMARK   8 BONDTYPE") == 0) {
				Model &model = mols.last().last().last();
				vector<string> vs = help::ssplit(line.substr(20), " ");
				Atom &a1 = *atom_number_to_atom[&model][stoi(vs[2])];
				Atom &a2 = *atom_number_to_atom[&model][stoi(vs[3])];
				auto &bond = a1.connect(a2);
				bond.set_bo(stoi(vs[0]));
				bond.set_bond_gaff_type(vs[1]);
			}
			else if (line.compare(0, 6, "CONECT") == 0) {
				const string ln = boost::algorithm::trim_right_copy(line.substr(6));
				dbgmsg("--" << ln << "--");
				vector<int> anum;
				for (int i = 0; i < ln.size(); i+=5) {
					const int atom_number = stoi(ln.substr(i, 5));
					dbgmsg(atom_number);
					anum.push_back(atom_number);
				}
				Molecule &molecule = mols.last(); // the CONECT words are valid for the whole molecule
				for (auto &assembly : molecule) {
					for (auto &model : assembly) {
						vector<int>::iterator it = anum.begin();
						auto atom_iterator1 = atom_number_to_atom[&model].find(*anum.begin());
						dbgmsg((atom_iterator1 == atom_number_to_atom[&model].end() 
							? "warning: atom1 not found (check CONECT words)" : ""));
						if (atom_iterator1 != atom_number_to_atom[&model].end()) {
							dbgmsg("found atom = " << boolalpha 
								<< (atom_iterator1 != atom_number_to_atom[&model].end()));
							Atom &a1 = *atom_iterator1->second;
							it++;
							while (it != anum.end()) {
								auto atom_iterator2 = atom_number_to_atom[&model].find(*it);
								dbgmsg((atom_iterator2 == atom_number_to_atom[&model].end() 
									? "warning: atom2 not found (check CONECT words)" : ""));
								if (atom_iterator2 != atom_number_to_atom[&model].end()) {
									dbgmsg("found atom = " << boolalpha 
										<< (atom_iterator2 != atom_number_to_atom[&model].end()));
									Atom &a2 = *atom_iterator2->second;
									dbgmsg(a1 << " " << a2);
									a1.connect(a2);
								}
								it++;
							}
						}
					}
				}
				connect_bonds(get_bonds_in(molecule.get_atoms()));
			}
		}
	}
};
