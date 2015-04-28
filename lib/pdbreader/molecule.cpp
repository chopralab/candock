#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "molecule.hpp"
#include "atomtype.hpp"
#include "bondtype.hpp"
#include "bond.hpp"
#include "geom3d/geom3d.hpp"
#include "fragmenter/fragmenter.hpp"
#include "fragmenter/unique.hpp"
#include "openmm/forcefield.hpp"
using namespace std;

namespace Molib {
	MolGraph create_graph(const AtomVec &atoms) {
		return MolGraph(atoms, true);
	}
	MolGraph create_graph(const AtomSet &atoms) {
		return MolGraph(atoms, true);
	}

	BondSet get_bonds_in(const AtomSet &atoms, bool in) {
		BondSet bonds;
		for (auto &patom1 : atoms) {
			for (auto &pbond : patom1->get_bonds()) {
				Atom &atom2 = pbond->second_atom(*patom1);
				if ((in || !atoms.count(&atom2)) 
					&& (!in || atoms.count(&atom2))) {
					bonds.insert(pbond);
				}
			}
		}
		return bonds;
	}
	BondSet get_bonds_in(const AtomVec &atoms, bool in) {
		AtomSet a(atoms.begin(), atoms.end());
		return get_bonds_in(a, in);
	}
	set<int> get_idatm_types(const Molecules &mols, set<int> previous) {
		for (auto &molecule : mols)
		for (auto &pa : molecule.get_atoms()) {
			previous.insert(pa->idatm_type());
		}
		return previous;
	}
	void Molecule::undo_mm_specific() {
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model) {
			for (auto &residue : chain) {
				if (residue.resn() == "HSD") residue.set_resn("HIS");
				if (residue.resn() == "CYX") residue.set_resn("CYS");
				if (residue.resn().size() == 4) residue.set_resn(residue.resn().substr(1)); // NALA -> ALA
			}
		}
	}

	Atom* get_closest_atom_of(const Atom &atom1, const AtomVec &neighbors, const string &atom_name) {
		double min_dist = HUGE_VAL;
		Atom *patom2 = nullptr;
		for (auto &pa : neighbors) {
			//~ if (!atom1.is_adjacent(*pa) && pa->atom_name() == atom_name) {
			if (pa->atom_name() == atom_name) {
				double d = atom1.crd().distance(pa->crd());
				if (d < min_dist) {
					min_dist = d;
					patom2 = pa;
				}
			}
		}
		return patom2;
	}
	void Molecule::prepare_for_mm(const OMMIface::ForceField &ff, MolGrid &grid) {
		/* Rename some residues
		 */
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if (residue.resn() == "HIS") residue.set_resn("HID");
			for (auto &atom : residue) {
				if (atom.atom_name() == "OXT") { // not foolproof, sometimes OXT is missing !!!
					residue.set_resn("C" + residue.resn()); // last residue is renamed to CALA, CGLY,...
					break;
				}
			}
		}
		/* Add bonds inside residues to protein atoms according to the topology file.
		 */
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model) {
			for (auto &residue : chain) {
				if (!ff.residue_topology.count(residue.resn())) 
					throw Error("die : cannot find topology for residue " + residue.resn());
				dbgmsg("residue topology for residue " << residue.resn());
				const OMMIface::ForceField::ResidueTopology &rtop = ff.residue_topology.at(residue.resn());
				for (auto &atom1 : residue) {
					if(rtop.bond.count(atom1.atom_name())) {
						for (auto &atom2 : residue) {
							if (rtop.bond.at(atom1.atom_name()).count(atom2.atom_name())) {
								atom1.connect(atom2);
							}
						}
					}
				}
			}
		}
		/* Add main chain and disulfide bonds
		 */
		for (auto &patom1 : this->get_atoms()) {
			auto &atom1 = *patom1;
			Atom *patom2 = nullptr;
#ifndef NDEBUG
			if (atom1.atom_name() == "N" || atom1.atom_name() == "C") {
				dbgmsg("neighbors of atom " << atom1.atom_name() << " " 
					<< atom1.atom_number() << " " << atom1.crd() << " are : ");
				for (auto &patom2 : grid.get_neighbors(atom1, 1.6)) {
					auto &atom2 = *patom2;
					dbgmsg("    neighbor : " << atom2.atom_name() << " " 
						<< atom2.atom_number() << " " << atom2.crd());
				}
			}
#endif
			if (atom1.atom_name() == "SG") {
				patom2 = get_closest_atom_of(atom1, grid.get_neighbors(atom1, 2.5), "SG");
			} else if (atom1.atom_name() == "N") {
				patom2 = get_closest_atom_of(atom1, grid.get_neighbors(atom1, 1.4), "C");
			}
			if (patom2) {
				auto &atom2 = *patom2;
				Residue &residue1 = const_cast<Residue&>(atom1.br());
				Residue &residue2 = const_cast<Residue&>(atom2.br());
				if (&residue1 != &residue2) {
					atom1.connect(atom2);
					dbgmsg("added inter-residue bond between atom " << atom1.atom_name() << " " 
						<< atom1.atom_number() << " and atom " << atom2.atom_name() << " " 
						<< atom2.atom_number());
					if (atom1.atom_name() == "SG") {
						residue1.set_resn("CYX");
						residue2.set_resn("CYX");
					}
				}
			}
		}
		/* Find terminal residues
		 */
		for (auto &patom : this->get_atoms()) {
			auto &atom = *patom;
			Residue &residue = const_cast<Residue&>(atom.br());
			if (residue.resn().size() == 3) {
				if (atom.atom_name() == "N" && !atom.is_adjacent("C")) {
					residue.set_resn("N" + residue.resn()); // first residue is renamed to NALA, NGLY,...
				} else if (atom.atom_name() == "C" && !atom.is_adjacent("N")) {
						residue.set_resn("C" + residue.resn()); // IF NOT ALREADY, rename last residue to CALA, CGLY,...
				}
			}
		}
		dbgmsg("MOLECULE AFTER PREPARING FOR MOLECULAR MECHANICS" << endl << *this);
	}
	AtomVec NRset::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		AtomVec atoms;
		for (auto &mols : *this) {
			auto ret = mols.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}
	AtomVec Molecules::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		AtomVec atoms;
		for (auto &molecule : *this) {
			auto ret = molecule.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}
	AtomVec Molecule::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		AtomVec atoms;
		for (auto &assembly : *this) {
			auto ret = assembly.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}
	AtomVec Assembly::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		AtomVec atoms;
		for (auto &model : *this) {
			if (model_number == -1 || model.number() == model_number) {
				auto ret = model.get_atoms(chain_ids, rest);
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}
	AtomVec Model::get_atoms(const string &chain_ids, const Residue::res_type &rest) const {
		AtomVec atoms;
		for (auto &chain : *this) {
			if (chain_ids == "" || chain_ids.find(chain.chain_id()) != string::npos) {
				auto ret = chain.get_atoms(rest);
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}
	AtomVec Chain::get_atoms(const Residue::res_type &rest) const {
		AtomVec atoms;
		for (auto &residue : *this) {
			if (rest == Residue::res_type::notassigned || residue.rest() == rest) {
				auto ret = residue.get_atoms();
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}
	AtomVec Residue::get_atoms() const {
		AtomVec atoms;
		for (auto &atom : *this) {
			atoms.push_back(&atom);
		}
		return atoms;
	}
	double compute_rmsd(const MolGraph &g1, const MolGraph &g2, 
		const MolGraph::Matches &m) {
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
					Atom &a1 = g1[vertices1[i]];
					Atom &a2 = g2[vertices2[i]];
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
	double Molecule::compute_rmsd(const Molecule &molecule) const {
		dbgmsg("calculate rmsd between two conformations of the same \
			molecule (can do symmetric molecules such as benzene, etc.)");
		MolGraph g1 = create_graph(this->get_atoms());
		dbgmsg("g1 = " << endl << g1);
		MolGraph g2 = create_graph(molecule.get_atoms());
		dbgmsg("g2 = " << endl << g2);
		if (g1.size() != g2.size()) 
			throw Error("die : RMSD can only be calculated for two conformations of the same molecule");
		MolGraph::Matches m = g1.match(g2);
		return Molib::compute_rmsd(g1, g2, m);
	}
	double Molecule::max_dist() const {
		dbgmsg("calculate max distance between atoms in an AtomSet ...");
		double md_sq = 0.0;
		AtomVec atoms = this->get_atoms();
		for (int i = 0; i < atoms.size(); ++i)
			for (int j = i + 1; j < atoms.size(); ++j) {
				const double dist_sq = atoms[i]->crd().distance_sq(atoms[j]->crd());
				if (dist_sq > md_sq) md_sq = dist_sq;
			}
		dbgmsg("max_dist between atoms of atoms in AtomSet = " << sqrt(md_sq));
		return sqrt(md_sq);
	}

	Molecules& Molecules::compute_seeds(const string &seeds_file) { 
		Unique u(seeds_file); 
		for (auto &molecule : *this) 
			molecule.compute_seeds(u); 
		return *this; 
	}

	Molecule& Molecule::compute_rotatable_bonds() {
		dbgmsg("starting compute_rotatable_bonds");
		Fragmenter f(this->get_atoms());
		f.substitute_bonds(help::rotatable);
		// rotatable bond inside a ring is changed back to non-rotatable
		Rings rings = f.identify_fused_rings();
		//~ Rings rings = f.identify_rings();
		for (auto &ring : rings) {
			for (auto &pbond : get_bonds_in(ring)) {
				pbond->set_rotatable("");
			}
		}
		dbgmsg("MOLECULE AFTER COMPUTING ROTATABLE BONDS" << endl << *this);
		return *this;
	}

	Molecule& Molecule::compute_overlapping_rigid_segments() {
		for (auto &assembly : *this)
		for (auto &model : assembly) {
			Fragmenter f(model.get_atoms());
			model.set_rigid(f.identify_overlapping_rigid_segments(model.get_atoms()));
		}
		dbgmsg("MOLECULE AFTER COMPUTING OVERLAPPING RIGID SEGMENTS" 
			<< endl << *this);
		return *this;
	}

	Molecule& Molecule::compute_seeds(Unique &u) {
		for (auto &assembly : *this)
		for (auto &model : assembly) {
			Fragmenter f(model.get_atoms());
			model.set_seeds(f.identify_seeds(model.get_rigid(), u));
		}
		return *this;
	}

	Molecule& Molecule::regenerate_bonds(const Molecule &template_molecule) {
		// copied molecule must be a fragment of template molecule...
		auto &molecule = *this;
		Molecule::const_iterator itt = template_molecule.begin();
		for (Molecule::iterator it = molecule.begin(); 
			it != molecule.end(); ++it, ++itt) {
			auto &assembly = *it;
			auto &template_assembly = *itt;
			for (auto it2 = assembly.begin(), itt2 = template_assembly.begin();
				it2 != assembly.end(); ++it2, ++itt2) {
				auto &model = *it2;
				auto &template_model = *itt2;
				model.regenerate_bonds(template_model);
			}
		}
		return *this;
	}

	Model& Model::regenerate_bonds(const Model &template_model) {
		// copied molecule must be a fragment of template molecule...
		map<int, Atom*> atom_number_to_atom;
		for (auto &patom : template_model.get_atoms()) {
			atom_number_to_atom[patom->atom_number()] = patom;
		}
		map<const Atom*, Atom*> atom1_to_copy1;
		for (auto &patom : this->get_atoms()) {
			patom->clear(); // clear bondees as they refer to the template molecule
			patom->clear_bonds(); // clear bonds : fixes double CONECT words bug in PDB output
			atom1_to_copy1[atom_number_to_atom.at(patom->atom_number())] = patom;
		}
		for (auto &kv : atom1_to_copy1) { // regenerate bonds
			const Atom &atom1 = *kv.first;
			Atom &copy1 = *kv.second;
			for (auto &atom2 : atom1) {
				if (atom1_to_copy1.count(&atom2)) {
					// copying of bonds does not preserve bond properties !!!
					Atom &copy2 = *atom1_to_copy1.at(&atom2);
					copy1.connect(copy2);
				}
			}
		}
		connect_bonds(get_bonds_in(this->get_atoms()));
		return *this;
	}

	Molecule& Molecule::filter(const unsigned int hm, const string &chain_ids) {
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model)
			if (chain_ids == "" || chain_ids.find(chain.chain_id()) != string::npos)
				for (auto &residue : chain)
					if ((hm & Residue::protein) && residue.rest() == Residue::protein
						|| (hm & Residue::nucleic) && residue.rest() == Residue::nucleic
						|| (hm & Residue::ion) && residue.rest() == Residue::ion
						|| (hm & Residue::water) && residue.rest() == Residue::water
						|| (hm & Residue::hetero) && residue.rest() == Residue::hetero) {
					} else
						chain.erase(Residue::res_pair(residue.resi(), residue.ins_code()));
			else
				model.erase(chain.chain_id());
	}
	
	Molecule& Molecule::compute_idatm_type() { 
		AtomType::compute_idatm_type(*this);
		return *this;
	}
	Molecule& Molecule::refine_idatm_type() { 
		AtomType::refine_idatm_type(*this);
		return *this;
	}
	Molecule& Molecule::compute_gaff_type() {
		AtomType::compute_gaff_type(*this);
		return *this;
	}
	Molecule& Molecule::compute_ring_type() {
		AtomType::compute_ring_type(*this);
		return *this;
	}
	Molecule& Molecule::compute_bond_gaff_type() {
		BondOrder::compute_bond_gaff_type(*this);
		return *this;
	}
	Molecule& Molecule::compute_bond_order() {
		BondOrder::compute_bond_order(*this);
		return *this;
	}
	int Atom::get_num_hydrogens() const {
		const Atom &atom = *this;
		int num_h = 0;
		for (auto &atom2 : atom)
			if (atom2.element() == Element::H) 
				num_h++;
		return num_h;
	}
	Bond& Atom::connect(Atom &a2) {
		Atom &a1 = *this;
		if (!a1.is_adjacent(a2)) {
			dbgmsg("connecting atoms " << a1.atom_number() << " and " << a2.atom_number());
			a1.add(&a2);
			a2.add(&a1);
			a1.insert_bond(a2, new Bond(&a1, &a2)); // insert if not exists
			return *a2.insert_bond(a1, a1.get_shared_ptr_bond(a2)); // insert if not exists
		}
		return a1.get_bond(a2); // if already connected return existing bond
	}
	void Atom::set_members(const string &str) {
		// set atomic penalty scores
		boost::smatch m;
		dbgmsg("set members for atom " << this->atom_name() << " " << this->atom_number());
		if (boost::regex_search(str, m, boost::regex("aps=(\\S+)"))) {
			if (m[1].matched) {
				string aps_str = m[1].str();
				dbgmsg("aps_str = " << aps_str);
				boost::match_results<string::const_iterator> matches;
				string::const_iterator begin = aps_str.begin(), end = aps_str.end();
				while (boost::regex_search(begin, end, matches, boost::regex("\\{(\\d+,\\d+)\\}"))) {
					string s(matches[1].first, matches[1].second);
					auto vec = help::ssplit(s, ",");
					int val = stoi(vec[0]);
					int aps = stoi(vec[1]);
					__aps[val] = aps;
					begin = matches[1].second;
				}
			}
		}
		// set properties
		if (boost::regex_search(str, m, boost::regex("prop=(\\S+)"))) {
			if (m[1].matched) {
				string prop_str = m[1].str();
				dbgmsg("prop = " << prop_str);
				boost::match_results<string::const_iterator> matches;
				string::const_iterator begin = prop_str.begin(), end = prop_str.end();
				while (boost::regex_search(begin, end, matches, boost::regex("\\{(\\w+,\\d+)\\}"))) {
					string s(matches[1].first, matches[1].second);
					auto vec = help::ssplit(s, ",");
					string val = vec[0];
					int count = stoi(vec[1]);
					__smiles_prop[val] = count;
					begin = matches[1].second;
				}
			}
		}
		// set gaff type
		if (boost::regex_search(str, m, boost::regex("gaff=([^,$]+)"))) {
			if (m[1].matched) {
				dbgmsg(m[1].str());
				set_gaff_type(m[1].str());
			}
		}
		// set idatm type
		if (boost::regex_search(str, m, boost::regex("idatm=([^,$]+)"))) {
			if (m[1].matched) {
				dbgmsg(m[1].str());
				set_idatm_type(m[1].str());
			}
		}
	}

	bool Atom::compatible(const Atom &other) const {
		/* this is "smiles" atom and other is real atom */
		if (!boost::regex_search(other.get_label(), boost::regex(this->get_label())))
			return false;
		for (auto &kv : this->__smiles_prop) {
			auto &p = kv.first;
			auto &cnt = kv.second;
			int other_cnt = other.compute_num_property(p);
			if (!((cnt == -1 && other_cnt > 0) || cnt == other_cnt))
				return false;
		}
		return true;
	}
	int Atom::compute_num_property(const string &prop) const {
		if (prop == "B") return this->size();
		else if (prop == "H") return this->get_num_hydrogens();
		else if (prop == "sb" || prop == "db" || prop == "tb" ||
			//~ prop == "SB" || prop == "DB" || prop == "TB" || prop == "AB") 
			prop == "SB" || prop == "DB" || prop == "DL") 
			return this->get_num_bond_with_bond_gaff_type(prop);
		else if (prop.substr(0,2) == "AR" || prop.substr(0,2) == "RG"
			 || prop.substr(0,2) == "NG" || prop.substr(0,2) == "ag")
			return this->get_num_property(prop);
	}
	int Atom::get_num_bond_with_bond_gaff_type(const string &prop) const {
		const Atom &atom = *this;
		int num_bwp = 0;
		for (auto &pbond : atom.get_bonds())
			if (pbond->get_bond_gaff_type() == prop) 
				num_bwp++;
		return num_bwp;
	}

	Molecule& Molecule::compute_hydrogen() {
		Molecule &molecule = *this;
		try {
			int max_atom_number = 0;
			for (auto &patom : molecule.get_atoms()) 
				if (patom->atom_number() > max_atom_number)
					max_atom_number = patom->atom_number();
			for (auto &patom : molecule.get_atoms()) {
				Atom &atom = *patom;
				Residue &residue = const_cast<Residue&> (atom.br());
				if (atom.element() != Element::H) { // don't visit "just added" hydrogens
					int con = help::infoMap.at(atom.idatm_type_unmask()).substituents;
			
					/* EXCEPTION for 3-substituted P */
					if (atom.idatm_type_unmask() == "P" && atom.size() == 3) 
						con = 3;
					/* ***************************** */

					int num_h = con - atom.size();
					dbgmsg("computing hydrogens for " << atom.idatm_type_unmask() << " con = "
						<< con << " atom.size() = " << atom.size());
					if (num_h > 0) {
						dbgmsg("computing missing hydrogens for molecule " << molecule.name()
							<< " atom " << atom << " con = " << con
							<< " atom.size() = " << atom.size() << " number of added hydrogens = "
							<< num_h);
						// add dummy hydrogen atoms with zero coordinates
						for (int i = 0; i < num_h; ++i) {
							const string idatm_type = (atom.element() == Element::C 
								? "HC" : "H");
							dbgmsg("idatm_type = " << idatm_type);
							dbgmsg("idatm_mask = " << help::idatm_mask.at(idatm_type));
							Atom &hatom = residue.add(new Atom(++max_atom_number, "H", 
								Geom3D::Coordinate(), help::idatm_mask.at(idatm_type)));
							atom.connect(hatom);
							dbgmsg("added hydrogen");
						}
					} else if (num_h < 0) {
						int h_excess = abs(num_h);
						dbgmsg("deleting excess hydrogens for molecule " << molecule.name()
							<< " because according to IDATM type : " 
							<< atom.idatm_type_unmask() << " this atom : "
							<< atom << " should have " << h_excess 
							<< " less hydrogens!");
						// deleting hydrogens
						for (int i = 0; i < atom.size(); ++i) {
							auto &bondee = atom[i];
							if (bondee.element() == Element::H) {
								if (h_excess-- > 0) {
									atom.erase(i--);
									auto &shpbond = atom.get_shared_ptr_bond(bondee);
									const Bond &deleted_bond = *shpbond;
									Bond::erase_stale_refs(deleted_bond, atom.get_bonds()); // delete references 
									dbgmsg("shared_count1 = " << shpbond.use_count());
									atom.erase_bond(bondee);
									dbgmsg("shared_count2 = " << shpbond.use_count());
									bondee.erase_bond(atom);
									dbgmsg("shared_count3 = " << shpbond.use_count());
									dbgmsg("residue before erasing hydrogen " << bondee.atom_number() << " for molecule " 
										<< molecule.name() << endl << residue);
									residue.erase(bondee.atom_number());
									dbgmsg("residue after erasing hydrogen " << bondee.atom_number() << " for molecule " 
										<< molecule.name() << endl << residue);
								}
							}
						}
						if (h_excess > 0)
							throw Error("die : deleting of excess hydrogens failed for molecule "
								+ molecule.name());
					}
				}
			}
			// connect the new hydrogen bonds with the rest of the bond-graph
			erase_bonds(get_bonds_in(molecule.get_atoms()));
			connect_bonds(get_bonds_in(molecule.get_atoms()));
			dbgmsg("MOLECULE AFTER COMPUTING HYDROGENS " << endl << *this);
			dbgmsg("BONDS AFTER COMPUTING HYDROGENS " << endl << get_bonds_in(this->get_atoms()));
		} catch (exception& e) {
			cerr << "errmesg : " << e.what() << " for molecule = " << molecule.name() << endl;
		}
		return *this;
	}
	Molecule& Molecule::erase_hydrogen() {
		Molecule &molecule = *this;
		for (auto &patom : molecule.get_atoms()) {
			Atom &atom = *patom;
			Residue &residue = const_cast<Residue&> (atom.br());
			if (atom.element() != Element::H) { // don't visit "just erased" hydrogens
				dbgmsg("deleting all hydrogens for molecule " << molecule.name()
					<< " for this atom : " << atom);
				// deleting hydrogens
				for (int i = 0; i < atom.size(); ++i) {
					auto &bondee = atom[i];
					if (bondee.element() == Element::H) {
						atom.erase(i--);
						atom.erase_bond(bondee);
						residue.erase(bondee.atom_number());
					}
				}
			}
		}
		dbgmsg("MOLECULE AFTER ERASING HYDROGENS " << endl << *this);
		return *this;
	}
	

	ostream& operator<< (ostream& stream, const Atom& a) {
		if (!a.__br) {
			stream << "ATOM " << a.get_label() << " " << a.atom_number() << endl;
			return stream;
		}
		stream << setw(6) << left << (a.br().rest() == Molib::Residue::hetero 
			|| a.br().rest() == Molib::Residue::ion 
			|| a.br().rest() == Molib::Residue::water ? "HETATM" : "ATOM");
		stream << setw(5) << right << a.atom_number();
		stream << setw(1) << " ";
		stream << setw(4) << left << (a.atom_name().size() < 4 ? " " + a.atom_name() : a.atom_name());
		stream << setw(1) << " ";
		stream << setw(3) << right << a.br().resn();
		stream << setw(1) << " ";
		stream << setw(1) << a.br().br().chain_id();
		stream << setw(4) << right << a.br().resi();
		stream << setw(1) << a.br().ins_code();
		stream << setw(27) << a.crd().pdb();
		stream << setw(6) << setprecision(2) << fixed << right << 1.0;
		stream << setw(6) << setprecision(2) << fixed << right << 1.0;
		stream << setw(12) << right << a.element();
		stream << setw(2) << "  ";
		stream << setw(5) << right << a.idatm_type_unmask();
		stream << setw(5) << right << a.gaff_type();
		stringstream ss;
		ss << " aps={";
		for (auto it = a.__aps.begin(); it != a.__aps.end(); ++it) {
			auto it_plus_one = it;
			ss << "{" << it->first << "," << it->second << "}" 
				<< (++it_plus_one == a.__aps.end() ? "" : ",");
		}
		stream << ss.str() << "}"; 
		ss.str("");
		ss << " prop={";
		for (auto it = a.__smiles_prop.begin(); it != a.__smiles_prop.end(); ++it) {
			auto it_plus_one = it;
			ss << "{" << it->first << "," << it->second << "}" 
				<< (++it_plus_one == a.__smiles_prop.end() ? "" : ",");
		}
		stream << ss.str() << "}"; 
		stream << endl;
		return stream;
	}
	ostream& operator<< (ostream& stream, const Residue& r) {
		for (auto &atom : r) { 
			stream << atom;
		} 
		return stream;
	}
	ostream& operator<< (ostream& stream, const Chain& c) {
		for (auto &residue : c) { 
			stream << residue;
		}
		stream << "TER" << endl;
		return stream;
	}

	ostream& operator<< (ostream& stream, const Model& m) {
		stream << setw(6) << left << "MODEL" << setw(4) << " " << setw(4) 
			<< right << m.number() << endl;
		for (auto &r1 : m.__remarks) {
			const int &remark_number = r1.first.first;
			const Residue::res_tuple2 &ligand = r1.first.second;
			stream << setw(6) << left << "REMARK" << setw(4) << right << remark_number 
				<< " " << std::get<0>(ligand) << ":" << std::get<1>(ligand) 
				<< ":" << std::get<2>(ligand) << ":" << std::get<3>(ligand);
			for (auto &rpair : r1.second) {
				stream <<  "," << std::get<0>(rpair.first) << ":" 
					<< std::get<1>(rpair.first) << ":" 
					<< std::get<2>(rpair.first) << ":" 
					<< std::get<3>(rpair.first);
				stream <<  "=" << std::get<0>(rpair.second) << ":" 
					<< std::get<1>(rpair.second) << ":" 
					<< std::get<2>(rpair.second) << ":" 
					<< std::get<3>(rpair.second);
			}
			stream << endl;
		}
		for (auto &chain : m) {
			stream << chain;
		}
		for (auto &kv : m.get_rigid()) {
			const string &nm = kv.first;
			const set<AtomSet> &same_name = kv.second;
			for (auto &rigid : same_name) {
				stream << "REMARK   8 RIGID " << nm;
				for (auto &atom : rigid) stream << " " << atom->atom_number();
				stream << endl;
			}
		}
		for (auto &kv : m.get_seeds()) {
			const string &nm = kv.first;
			const set<AtomSet> &same_name = kv.second;
			for (auto &seed : same_name) {
				stream << "REMARK   8 SEED " << nm;
				for (auto &atom : seed) stream << " " << atom->atom_number();
				stream << endl;
			}
		}
		for (auto &pbond : get_bonds_in(m.get_atoms())) {
			Bond &bond = *pbond;
			if (bond.is_rotatable()) {
				stream << "REMARK   8 ROTA " << bond.get_rotatable() 
					<< " " << bond.atom1().atom_number() 
					<< " " << bond.atom2().atom_number() 
					<< endl;
			}
		}
		for (auto &pbond : get_bonds_in(m.get_atoms())) {
			Bond &bond = *pbond;
			stream << "REMARK   8 BONDTYPE " << bond.get_bo() 
				<< " " << bond.get_bond_gaff_type() 
				<< " " << bond.atom1().atom_number() 
				<< " " << bond.atom2().atom_number() 
				<< endl;
		}
		for(auto &chain : m)
		for(auto &residue : chain) {
			// don't write conect for standard residues
			if (!help::standard_residues.count(residue.resn())) {
				for(auto &atom : residue) {
					for (auto &adj_a : atom) {
						int bond_order = atom.get_bond(adj_a).get_bo();
						if (bond_order == 0) bond_order = 1;
						for (int bo = 0; bo < bond_order; ++bo) {
							stream << "CONECT" << setw(5) << right 
								<< atom.atom_number() << setw(5) 
								<< right << adj_a.atom_number() << endl;
						}
					}
				}
			}
		}
		stream << "ENDMDL" << endl;
		return stream;
	}
	ostream& operator<< (ostream& stream, const Assembly& a) {
		stream << setw(6) << left << "REMARK   6 " << a.name() << " " << a.number() << endl;
		for (auto &model : a) { 
			stream << model;
		}
		stream << "REMARK   6 END" << endl;
		return stream;
	}
	ostream& operator<< (ostream& stream, const Molecule& m) {
		stream << setw(6) << left << "REMARK   5 MOLECULE " << m.name() << endl;
		for (auto &assembly : m) {
			stream << assembly;
		}
		stream << setw(6) << left << "REMARK   5 END" << endl;
		return stream;
	}
	ostream& operator<< (ostream& stream, const Molecules& m) {
		stream << setw(6) << left << "REMARK   4 NRPDB " << m.name() << endl;
		for (auto &molecule : m) {
			stream << molecule;
		}
		stream << setw(6) << left << "REMARK   4 END" << endl;
		return stream;
	}
	ostream& operator<< (ostream& stream, const NRset& m) {
		for (auto &molecules : m) {
			stream << molecules;
		}
		return stream;
	}

	ostream& operator<< (ostream& stream, const AtomSet&atoms) {
		for (auto &patom : atoms) stream << *patom << endl;
		return stream;
	}
	ostream& operator<< (ostream& stream, const AtomVec&atoms) {
		for (auto &patom : atoms) stream << *patom << endl;
		return stream;
	}

};
