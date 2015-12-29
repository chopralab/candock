#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
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
#include "modeler/forcefield.hpp"
using namespace std;

namespace Molib {
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
				if (residue.resn() == "HIP") residue.set_resn("HIS");
				if (residue.resn() == "CYX") residue.set_resn("CYS");
				if (residue.resn().size() == 4) residue.set_resn(residue.resn().substr(1)); // NALA -> ALA
			}
		}
	}

	Atom* get_closest_atom_of(const Atom &atom1, const Atom::Vec &neighbors, const string &atom_name) {
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

	void Molecule::prepare_for_mm(const OMMIface::ForceField &ff, Atom::Grid &grid) {
		/* Rename some residues
		 */
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if (residue.resn() == "HIS") residue.set_resn("HIP");
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
				for (auto &patom2 : grid.get_neighbors(atom1.crd(), 1.6)) {
					auto &atom2 = *patom2;
					dbgmsg("    neighbor : " << atom2.atom_name() << " " 
						<< atom2.atom_number() << " " << atom2.crd());
				}
			}
#endif
			if (atom1.atom_name() == "SG") {
				patom2 = get_closest_atom_of(atom1, grid.get_neighbors(atom1.crd(), 2.5), "SG");
			} else if (atom1.atom_name() == "N") {
				patom2 = get_closest_atom_of(atom1, grid.get_neighbors(atom1.crd(), 1.4), "C");
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
	
	
	Molecule::Vec Molecules::get_molecules(const Residue::res_type &rest) const {
		Molecule::Vec molecules;
		for (auto &molecule : *this) {
			Residue &first = molecule.first().first().first().first();
			if (first.rest() == rest)
				molecules.push_back(&molecule);
		}
		return molecules;
	}

	void NRset::jiggle() {
		srand((unsigned)time(NULL));
		for (auto &patom : this->get_atoms()) {
			dbgmsg("before jiggle crd = " << patom->crd());
			Geom3D::Coordinate dcrd(((double)rand() / (double)RAND_MAX) / 1000, 
				((double)rand() / (double)RAND_MAX) / 1000,
				((double)rand() / (double)RAND_MAX) / 1000);
			patom->set_crd(patom->crd() + dcrd);
			dbgmsg("after jiggle crd = " << patom->crd());
		}
	}

	Molecule::Vec NRset::get_molecules(const Residue::res_type &rest) const {
		Molecule::Vec molecules;
		for (auto &mols : *this) {
			Residue &first = mols.first().first().first().first().first();
			auto ret = mols.get_molecules(rest);
			molecules.insert(molecules.end(), ret.begin(), ret.end());
		}
		return molecules;
	}

	Atom::Vec NRset::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		Atom::Vec atoms;
		for (auto &mols : *this) {
			auto ret = mols.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}
	Atom::Vec Molecules::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		Atom::Vec atoms;
		for (auto &molecule : *this) {
			auto ret = molecule.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}
	Atom::Vec Molecule::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		Atom::Vec atoms;
		for (auto &assembly : *this) {
			auto ret = assembly.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}

	Geom3D::Point::Vec Molecule::get_crds(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const { 

		Geom3D::Point::Vec crds;
		for (auto &patom : this->get_atoms(chain_ids, rest, model_number)) 
			crds.push_back(patom->crd());
		return crds;
	}

	Geom3D::Point::Vec Molecules::get_crds(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const { 

		Geom3D::Point::Vec crds;
		for (auto &patom : this->get_atoms(chain_ids, rest, model_number)) 
			crds.push_back(patom->crd());
		return crds;
	}

	Atom::Vec Assembly::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		Atom::Vec atoms;
		for (auto &model : *this) {
			if (model_number == -1 || model.number() == model_number) {
				auto ret = model.get_atoms(chain_ids, rest);
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}

	Atom::Vec Model::get_atoms(const string &chain_ids, const Residue::res_type &rest) const {
		Atom::Vec atoms;
		for (auto &chain : *this) {
			if (chain_ids == "" || chain_ids.find(chain.chain_id()) != string::npos) {
				auto ret = chain.get_atoms(rest);
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}
	Atom::Vec Chain::get_atoms(const Residue::res_type &rest) const {
		Atom::Vec atoms;
		for (auto &residue : *this) {
			if (rest == Residue::res_type::notassigned || residue.rest() == rest) {
				auto ret = residue.get_atoms();
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}
	Atom::Vec Residue::get_atoms() const {
		Atom::Vec atoms;
		for (auto &atom : *this) {
			atoms.push_back(&atom);
		}
		return atoms;
	}

	double Molecule::compute_rmsd(const Molecule &molecule) const {
		dbgmsg("calculate rmsd between two conformations of the same \
			molecule (can do symmetric molecules such as benzene, etc.)");

		Atom::Graph g1 = Atom::create_graph(this->get_atoms());
		dbgmsg("g1 = " << endl << g1);

		Atom::Graph g2 = Atom::create_graph(molecule.get_atoms());
		dbgmsg("g2 = " << endl << g2);

		if (g1.size() != g2.size())
			throw Error("die : RMSD can only be calculated for two conformations of the same molecule");

		Atom::Graph::Matches m = g1.match(g2);
		
		// try calculating rmsd of each mapping of molecule to molecule ...
		double min_sum_squared = HUGE_VAL;

		for (auto &mv : m) {
			auto &vertices1 = mv.first;
			auto &vertices2 = mv.second;
			// match must stretch over the whole graph g1 (and g2)
			if (vertices2.size() != g1.size()) 
				throw Error("die : RMSD calculation aborted due to imperfect match"); 
			double sum_squared = 0;
			for (int i = 0; i < vertices2.size(); ++i) {
				Atom &a1 = g1[vertices1[i]];
				Atom &a2 = g2[vertices2[i]];
				sum_squared += a1.crd().distance_sq(a2.crd());
			}
			if (sum_squared < min_sum_squared)
				min_sum_squared = sum_squared;
		}

		return sqrt(min_sum_squared / g1.size());	
	}

	double Molecule::compute_rmsd_ord(const Molecule &other) const {
		return sqrt(Geom3D::compute_rmsd_sq(this->get_crds(), other.get_crds()));
	}

	Atom& Molecule::get_center_atom() const {
		Atom *center_atom = nullptr;
		Atom::Vec atoms = this->get_atoms();
		double min_d = HUGE_VAL;
		for (auto &patom1 : atoms) {
			double d = 0.0;
			for (auto &patom2 : atoms) {
				if (patom1 != patom2) {
					d += patom1->crd().distance(patom2->crd());
				}
			}
			if (d < min_d) {
				min_d = d;
				center_atom = patom1;
			}
		}
		if (!center_atom)
			throw Error("die : center atom cannot be found");
		return *center_atom;
	}

	double Molecule::max_dist() const {
		dbgmsg("calculate max distance between atoms in an Atom::Set ...");
		double md_sq = 0.0;
		Atom::Vec atoms = this->get_atoms();
		for (int i = 0; i < atoms.size(); ++i)
			for (int j = i + 1; j < atoms.size(); ++j) {
				const double dist_sq = atoms[i]->crd().distance_sq(atoms[j]->crd());
				if (dist_sq > md_sq) md_sq = dist_sq;
			}
		dbgmsg("max_dist between atoms of atoms in Atom::Set = " << sqrt(md_sq));
		return sqrt(md_sq);
	}
	double Molecule::max_dist(const Atom &atom) const {
		dbgmsg("calculate max distance between atom " << atom 
			<< " and other atoms in molecule ...");
		double md_sq = 0.0;
		Atom::Vec atoms = this->get_atoms();
		for (int i = 0; i < atoms.size(); ++i) {
			const double dist_sq = atom.crd().distance_sq(atoms[i]->crd());
			if (dist_sq > md_sq) md_sq = dist_sq;
		}
		return sqrt(md_sq);
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

	Molecules& Molecules::compute_overlapping_rigid_segments(const string &seeds_file) { 
		Unique u(seeds_file); 
		for (auto &molecule : *this) 
			molecule.compute_overlapping_rigid_segments(u); 
		return *this; 
	}

	Molecule& Molecule::compute_overlapping_rigid_segments(Unique &u) {
		for (auto &assembly : *this)
		for (auto &model : assembly) {
			Fragmenter f(model.get_atoms());
			model.set_rigid(f.identify_overlapping_rigid_segments(model.get_atoms(), u));
		}
		dbgmsg("MOLECULE AFTER COMPUTING OVERLAPPING RIGID SEGMENTS" 
			<< endl << *this);
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

		Molecule saved(*this);
		
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model)
			if (chain_ids == "" || chain_ids.find(chain.chain_id()) != string::npos) {
				for (auto &residue : chain)
					if ((hm & Residue::protein) && residue.rest() == Residue::protein
						|| (hm & Residue::nucleic) && residue.rest() == Residue::nucleic
						|| (hm & Residue::ion) && residue.rest() == Residue::ion
						|| (hm & Residue::water) && residue.rest() == Residue::water
						|| (hm & Residue::hetero) && residue.rest() == Residue::hetero) {
					} else {
						dbgmsg("erasing residue " << residue.resi() << " ins_code = " << residue.ins_code());
						chain.erase(Residue::res_pair(residue.resi(), residue.ins_code()));
					}
			}
			else {
				dbgmsg("erasing chain " << chain.chain_id());
				model.erase(chain.chain_id());
			}
			
		regenerate_bonds(saved);
		dbgmsg("out of filter");
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
				<< setw(6) << right << "ALRES" << " " << std::get<0>(ligand) << ":" << std::get<1>(ligand) 
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
		for (auto &rigid : m.get_rigid()) {
			stream << "REMARK   8 RIGID " << rigid.get_seed_id();
			for (auto &atom : rigid.get_core()) stream << " " << 'c' << atom->atom_number();
			for (auto &atom : rigid.get_join()) stream << " " << 'j' << atom->atom_number();
			stream << endl;
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

	/**
	 * Copy molecule with new coordinates
	 */
	Molecule::Molecule(const Molib::Molecule &rhs, const Geom3D::Point::Vec &crds) : Molecule(rhs) {

		auto atoms = get_atoms();
		for (int i = 0; i < atoms.size(); ++i) {
			atoms[i]->set_crd(crds[i]);
		}
	}

	string Molecule::print_complex(Molecule &ligand, Molecule &receptor, const double energy) {
		stringstream ss;
		ss << "REMARK   1 MINIMIZED COMPLEX OF " << ligand.name() << " AND " << receptor.name() << " WITH SCORE OF " << energy << endl;
		ss << "MODEL" << endl;
		int reenum = 0;
		for (auto &patom : receptor.get_atoms()) {
			ss << setw(66) << left << *patom;
			reenum = patom->atom_number();
		}
		ss << "TER" << endl;
		for (auto &patom : ligand.get_atoms()) {
			patom->set_atom_number(++reenum);
			ss << setw(66) << left << *patom;
		}
		//~ ss << get_bonds_in(ligand.get_atoms());
		ss << "ENDMDL" << endl;
		return ss.str();
	}


};
