#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "bond.hpp"
#include "hydrogens.hpp"
#include "geom3d/geom3d.hpp"
#include "fragmenter/fragmenter.hpp"
#include "fragmenter/unique.hpp"
#include "modeler/forcefield.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
#include "model.hpp"
#include "assembly.hpp"
#include "molecule.hpp"

using namespace std;

namespace Molib {
	
	string Molecule::get_chain_ids(const unsigned int hm) const {

		set<char> chain_ids;

		for (auto &presidue : this->get_residues()) {
			auto &residue = *presidue;
			if ((hm & Residue::protein) && residue.rest() == Residue::protein
				|| (hm & Residue::nucleic) && residue.rest() == Residue::nucleic
				|| (hm & Residue::ion) && residue.rest() == Residue::ion
				|| (hm & Residue::water) && residue.rest() == Residue::water
				|| (hm & Residue::hetero) && residue.rest() == Residue::hetero) {
				
				chain_ids.insert(residue.br().chain_id());
				
			}
		}
		if (chain_ids.empty()) throw Error("die : no protein in the receptor file?");
		string chains("");
		
		for (auto &ch : chain_ids) chains += ch;
		dbgmsg("these are protein chains in the receptor file = " << chains);
		return chains;
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

	void Molecule::rotate(const Geom3D::Matrix &rota, const bool inverse) {
		for (auto &assembly : *this) {
			assembly.rotate(rota, inverse);
		}
	}

	void Molecule::init_bio(bio_how_many bhm) {
		for(auto &i : __bio_rota) {
			int biomolecule_number = i.first;
			M0 &matrices = i.second;
			Assembly &last = add(new Assembly(biomolecule_number, "BIOLOGICAL ASSEMBLY"));
			last.init_bio(asym(), matrices, __bio_chain[biomolecule_number]);
			if (bhm == Molecule::first_bio) { // exit after initializing the first assembly
				break;
			}
		}
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
		dbgmsg("computing rmsd between molecule " << *this << " and " << other);
		return Geom3D::compute_rmsd(this->get_crds(), other.get_crds());
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
	
	Residue::Vec Molecule::get_residues() const {
		Residue::Vec residues;
		for (auto &assembly : *this)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			residues.push_back(&residue);
		}
		return residues;
	}



	void Molecule::prepare_for_mm(const OMMIface::ForceField &ffield, Atom::Grid &grid) {
		/* Rename some residues
		 */
		for (auto &presidue : this->get_residues()) {
			auto &residue = *presidue;
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
		for (auto &presidue : this->get_residues()) {
			auto &residue = *presidue;
			if (ffield.residue_topology.count(residue.resn())) {
				dbgmsg("residue topology for residue " << residue.resn());
				const OMMIface::ForceField::ResidueTopology &rtop = ffield.residue_topology.at(residue.resn());
				for (auto &atom1 : residue) {
					if(rtop.bond.count(atom1.atom_name())) {
						for (auto &atom2 : residue) {
							if (rtop.bond.at(atom1.atom_name()).count(atom2.atom_name())) {
								atom1.connect(atom2);
							}
						}
					}
				}
			} else {
				cerr << "note : cannot find topology for residue " << residue.resn() 
					<< " (this is normal if receptor file contains water, cofactors etc.)" << endl;
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
	
	ostream& operator<< (ostream& stream, const Molecule& m) {
		stream << setw(6) << left << "REMARK   5 MOLECULE " << m.name() << endl;
		for (auto &assembly : m) {
			stream << assembly;
		}
		stream << setw(6) << left << "REMARK   5 END" << endl;
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
