#include "moleculeinfo.hpp"
#include "forcefield.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "helper/error.hpp"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
using namespace std;

namespace OMMIface {
	ostream& operator<< (ostream& stream, const MoleculeInfo::Bonds &bonds) {
		for (auto &atom_pair : bonds) {
			auto &atom1 = *atom_pair.first;
			auto &atom2 = *atom_pair.second;
			stream << "bond between atoms " << atom1.atom_number() << " and " 
				<< atom2.atom_number() << endl;
		}
		return stream;
	}
	ostream& operator<< (ostream& stream, const MoleculeInfo::Angles &angles) {
		for (auto &angle : angles) {
			auto &atom1 = *get<0>(angle);
			auto &atom2 = *get<1>(angle);
			auto &atom3 = *get<2>(angle);
			stream << "angle between atoms " << atom1.atom_number() << " and " 
				<< atom2.atom_number() << " and " << atom3.atom_number() << endl;
		}
		return stream;
	}
	ostream& operator<< (ostream& stream, const MoleculeInfo::Dihedrals &dihedrals) {
		for (auto &dihedral : dihedrals) {
			auto &atom1 = *get<0>(dihedral);
			auto &atom2 = *get<1>(dihedral);
			auto &atom3 = *get<2>(dihedral);
			auto &atom4 = *get<3>(dihedral);
			stream << "dihedral between atoms " << atom1.atom_number() << " and " 
				<< atom2.atom_number() << " and " << atom3.atom_number() 
				<< " and " << atom4.atom_number() << endl;
		}
		return stream;
	}

	//~ MoleculeInfo& MoleculeInfo::get_molecule_info(const Molib::Molecule &molecule, 
					//~ const ForceField &ff) {
		//~ // residue topology
		//~ for (auto &assembly : molecule)
		//~ for (auto &model : assembly)
		//~ for (auto &chain : model) {
			//~ for (auto &residue : chain) {
				//~ if (!ff.residue_topology.count(residue.resn())) 
					//~ throw Error("die : cannot find topology for residue " + residue.resn());
				//~ dbgmsg("residue topology for residue " << residue.resn());
				//~ const ForceField::ResidueTopology &rtop = ff.residue_topology.at(residue.resn());
				//~ for (auto &atom : residue) {
					//~ if (!rtop.atom.count(atom.atom_name())) 
						//~ throw Error("die : cannot find atom " + atom.atom_name() 
								//~ + " in topology for residue " + residue.resn());
					//~ this->atom.push_back(&atom);
					//~ dbgmsg("crd = " << atom.crd() << " type = " << rtop.atom.at(atom.atom_name()));
					//~ this->atom_to_type.insert({&atom, rtop.atom.at(atom.atom_name())});
				//~ }
			//~ }
		//~ }
		//~ // set the bonds, angles, dihedrals
		//~ for (auto &assembly : molecule)
		//~ for (auto &model : assembly)
		//~ for (auto &chain : model)
		//~ for (auto &residue : chain)
		//~ for (auto &atom1 : residue) {
			//~ for (auto &atom2 : atom1) {
				//~ this->bond.push_back({&atom1, &atom2});
				//~ dbgmsg("info.bond = " << atom1.atom_number() << " " << atom2.atom_number());
				//~ for (auto &atom3 : atom2) {
					//~ if (&atom3 != &atom1) {
						//~ this->angle.push_back(make_tuple(&atom1, &atom2, &atom3));
						//~ dbgmsg("info.angle = " << atom1.atom_number() << " " << atom2.atom_number()
							//~ << " " << atom3.atom_number());
						//~ for (auto &atom4 : atom3) { // propers
							//~ if (&atom4 != &atom2 && &atom4 != &atom1) {
								//~ this->dihedral.push_back(make_tuple(&atom1, &atom2, &atom3, &atom4));
								//~ dbgmsg("info.dihedral (proper) = " << atom1.atom_number() << " " << atom2.atom_number()
									//~ << " " << atom3.atom_number() << " " << atom4.atom_number());
							//~ }
						//~ }
						//~ for (auto &atom4 : atom2) { // impropers
							//~ if (&atom4 != &atom3 && &atom4 != &atom1) {
								//~ this->dihedral.push_back(make_tuple(&atom3, &atom1, &atom2, &atom4));
								//~ dbgmsg("info.dihedral (improper) = " << atom3.atom_number() << " " << atom1.atom_number() 
									//~ << " " << atom2.atom_number() << " " << atom4.atom_number());
							//~ }
						//~ }
					//~ }
				//~ }
			//~ }
		//~ }
		//~ return *this;
	//~ }

	//~ MoleculeInfo& MoleculeInfo::get_molecule_info(const Molib::Molecule &molecule, 
					//~ const ForceField &ff) {
		//~ // residue topology
		//~ for (auto &assembly : molecule)
		//~ for (auto &model : assembly)
		//~ for (auto &chain : model) {
			//~ for (auto &residue : chain) {
				//~ if (!ff.residue_topology.count(residue.resn())) 
					//~ throw Error("die : cannot find topology for residue " + residue.resn());
				//~ dbgmsg("residue topology for residue " << residue.resn());
				//~ const ForceField::ResidueTopology &rtop = ff.residue_topology.at(residue.resn());
				//~ for (auto &atom : residue) {
					//~ if (!rtop.atom.count(atom.atom_name())) 
						//~ throw Error("die : cannot find atom " + atom.atom_name() 
								//~ + " in topology for residue " + residue.resn());
					//~ this->atom.push_back(&atom);
					//~ dbgmsg("crd = " << atom.crd() << " type = " << rtop.atom.at(atom.atom_name()));
					//~ this->atom_to_type.insert({&atom, rtop.atom.at(atom.atom_name())});
				//~ }
			//~ }
		//~ }
		//~ set<Molib::AtomPair> visited_bonds;
		//~ set<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*>> visited_angles;
		//~ set<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*, Molib::Atom*>> visited_dihedrals, visited_impropers;
		//~ // set the bonds, angles, dihedrals
		//~ for (auto &patom1 : molecule.get_atoms()) {
			//~ auto &atom1 = *patom1;
			//~ for (auto &atom2 : atom1) {
				//~ if (visited_bonds.count({&atom1, &atom2})) continue;
				//~ visited_bonds.insert({&atom2, &atom1});
				//~ 
				//~ this->bond.push_back({&atom1, &atom2});
				//~ dbgmsg("info.bond = " << atom1.atom_number() << " " << atom2.atom_number());
				//~ for (auto &atom3 : atom2) {
					//~ if (&atom3 != &atom1) {
						//~ if (visited_angles.count(make_tuple(&atom1, &atom2, &atom3))) continue;
						//~ visited_angles.insert(make_tuple(&atom3, &atom2, &atom1));
//~ 
						//~ this->angle.push_back(make_tuple(&atom1, &atom2, &atom3));
						//~ dbgmsg("info.angle = " << atom1.atom_number() << " " << atom2.atom_number()
							//~ << " " << atom3.atom_number());
						//~ for (auto &atom4 : atom3) { // propers
							//~ if (&atom4 != &atom2 && &atom4 != &atom1) {
								//~ if (visited_dihedrals.count(make_tuple(&atom1, &atom2, &atom3, &atom4))) continue;
								//~ visited_dihedrals.insert(make_tuple(&atom4, &atom3, &atom2, &atom1));
//~ 
								//~ this->dihedral.push_back(make_tuple(&atom1, &atom2, &atom3, &atom4));
								//~ dbgmsg("info.dihedral (proper) = " << atom1.atom_number() 
									//~ << " " << atom2.atom_number() << " " << atom3.atom_number() 
									//~ << " " << atom4.atom_number());
							//~ }
						//~ }
						//~ for (auto &atom4 : atom2) { // impropers
							//~ if (&atom4 != &atom3 && &atom4 != &atom1) {
								//~ if (visited_impropers.count(make_tuple(&atom3, &atom1, &atom2, &atom4)) 
									//~ || visited_impropers.count(make_tuple(&atom1, &atom3, &atom2, &atom4)) 
									//~ || visited_impropers.count(make_tuple(&atom1, &atom4, &atom2, &atom3)) 
									//~ || visited_impropers.count(make_tuple(&atom4, &atom1, &atom2, &atom3)) 
									//~ || visited_impropers.count(make_tuple(&atom4, &atom3, &atom2, &atom1)) 
									//~ || visited_impropers.count(make_tuple(&atom3, &atom4, &atom2, &atom1))) continue;
								//~ visited_impropers.insert(make_tuple(&atom4, &atom3, &atom2, &atom1));
//~ 
								//~ this->improper.push_back(make_tuple(&atom3, &atom1, &atom2, &atom4));
								//~ dbgmsg("info.improper (improper) (third atom is central) = " 
									//~ << atom3.atom_number() << " " << atom1.atom_number() 
									//~ << " " << atom2.atom_number() << " " << atom4.atom_number());
							//~ }
						//~ }
					//~ }
				//~ }
			//~ }
		//~ }
		//~ dbgmsg("BONDS ARE : " << endl << this->bond);
		//~ dbgmsg("ANGLES ARE : " << endl << this->angle);
		//~ dbgmsg("DIHEDRALS ARE : " << endl << this->dihedral);
		//~ dbgmsg("IMPROPERS ARE : " << endl << this->improper);
		//~ return *this;
	//~ }
	
	MoleculeInfo& MoleculeInfo::get_molecule_info(const Molib::Molecule &molecule, 
					const ForceField &ff) {
		// residue topology
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model) {
			for (auto &residue : chain) {
				if (!ff.residue_topology.count(residue.resn())) 
					throw Error("die : cannot find topology for residue " + residue.resn());
				dbgmsg("residue topology for residue " << residue.resn());
				const ForceField::ResidueTopology &rtop = ff.residue_topology.at(residue.resn());
				for (auto &atom : residue) {
					if (!rtop.atom.count(atom.atom_name())) 
						throw Error("die : cannot find atom " + atom.atom_name() 
								+ " in topology for residue " + residue.resn());
					this->atom.push_back(&atom);
					dbgmsg("crd = " << atom.crd() << " type = " << rtop.atom.at(atom.atom_name()));
					this->atom_to_type.insert({&atom, rtop.atom.at(atom.atom_name())});
				}
			}
		}
		set<Molib::AtomPair> visited_bonds;
		set<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*>> visited_angles;
		set<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*, Molib::Atom*>> visited_dihedrals, visited_impropers;
		// set the bonds, angles, dihedrals
		for (auto &patom1 : molecule.get_atoms()) {
			auto &atom1 = *patom1;
			for (auto &atom2 : atom1) {
				if (!visited_bonds.count({&atom1, &atom2})) {
					visited_bonds.insert({&atom2, &atom1});
					this->bond.push_back({&atom1, &atom2});
					dbgmsg("info.bond = " << atom1.atom_number() << " " << atom2.atom_number());
				}
				for (auto &atom3 : atom2) {
					if (&atom3 != &atom1) {
						if (!visited_angles.count(make_tuple(&atom1, &atom2, &atom3))) {
							visited_angles.insert(make_tuple(&atom3, &atom2, &atom1));
							this->angle.push_back(make_tuple(&atom1, &atom2, &atom3));
							dbgmsg("info.angle = " << atom1.atom_number() << " " << atom2.atom_number()
								<< " " << atom3.atom_number());
						}
						for (auto &atom4 : atom3) { // propers
							if (&atom4 != &atom2 && &atom4 != &atom1) {
								if (!visited_dihedrals.count(make_tuple(&atom1, &atom2, &atom3, &atom4))) {
									visited_dihedrals.insert(make_tuple(&atom4, &atom3, &atom2, &atom1));
									this->dihedral.push_back(make_tuple(&atom1, &atom2, &atom3, &atom4));
									dbgmsg("info.dihedral (proper) = " << atom1.atom_number() 
										<< " " << atom2.atom_number() << " " << atom3.atom_number() 
										<< " " << atom4.atom_number());
								}
							}
						}
						for (auto &atom4 : atom2) { // impropers
							if (&atom4 != &atom3 && &atom4 != &atom1) {
								if (!(visited_impropers.count(make_tuple(&atom3, &atom1, &atom2, &atom4)) 
									|| visited_impropers.count(make_tuple(&atom1, &atom3, &atom2, &atom4)) 
									|| visited_impropers.count(make_tuple(&atom1, &atom4, &atom2, &atom3)) 
									|| visited_impropers.count(make_tuple(&atom4, &atom1, &atom2, &atom3)) 
									|| visited_impropers.count(make_tuple(&atom4, &atom3, &atom2, &atom1)) 
									|| visited_impropers.count(make_tuple(&atom3, &atom4, &atom2, &atom1)))) {

									visited_impropers.insert(make_tuple(&atom4, &atom3, &atom2, &atom1));
									this->improper.push_back(make_tuple(&atom3, &atom1, &atom2, &atom4));
									dbgmsg("info.improper (improper) (third atom is central) = " 
										<< atom3.atom_number() << " " << atom1.atom_number() 
										<< " " << atom2.atom_number() << " " << atom4.atom_number());
								}
							}
						}
					}
				}
			}
		}
		//~ // set the bonds
		//~ for (auto &patom1 : molecule.get_atoms()) {
			//~ auto &atom1 = *patom1;
			//~ for (auto &atom2 : atom1) {
				//~ if (visited_bonds.count({&atom1, &atom2})) continue;
				//~ visited_bonds.insert({&atom2, &atom1});
				//~ 
				//~ this->bond.push_back({&atom1, &atom2});
				//~ dbgmsg("info.bond = " << atom1.atom_number() << " " << atom2.atom_number());
			//~ }
		//~ }
		//~ // set the angles
		//~ for (auto &patom1 : molecule.get_atoms()) {
			//~ auto &atom1 = *patom1;
			//~ for (auto &atom2 : atom1) {
				//~ for (auto &atom3 : atom2) {
					//~ if (&atom3 != &atom1) {
						//~ if (visited_angles.count(make_tuple(&atom1, &atom2, &atom3))) continue;
						//~ visited_angles.insert(make_tuple(&atom3, &atom2, &atom1));
//~ 
						//~ this->angle.push_back(make_tuple(&atom1, &atom2, &atom3));
						//~ dbgmsg("info.angle = " << atom1.atom_number() << " " << atom2.atom_number()
							//~ << " " << atom3.atom_number());
					//~ }
				//~ }
			//~ }
		//~ }
		//~ // set the dihedrals
		//~ for (auto &patom1 : molecule.get_atoms()) {
			//~ auto &atom1 = *patom1;
			//~ for (auto &atom2 : atom1) {
				//~ for (auto &atom3 : atom2) {
					//~ if (&atom3 != &atom1) {
						//~ for (auto &atom4 : atom3) { // propers
							//~ if (&atom4 != &atom2 && &atom4 != &atom1) {
								//~ if (visited_dihedrals.count(make_tuple(&atom1, &atom2, &atom3, &atom4))) continue;
								//~ visited_dihedrals.insert(make_tuple(&atom4, &atom3, &atom2, &atom1));
//~ 
								//~ this->dihedral.push_back(make_tuple(&atom1, &atom2, &atom3, &atom4));
								//~ dbgmsg("info.dihedral (proper) = " << atom1.atom_number() 
									//~ << " " << atom2.atom_number() << " " << atom3.atom_number() 
									//~ << " " << atom4.atom_number());
							//~ }
						//~ }
					//~ }
				//~ }
			//~ }
		//~ }
		//~ // set the impropers
		//~ for (auto &patom1 : molecule.get_atoms()) {
			//~ auto &atom1 = *patom1;
			//~ for (auto &atom2 : atom1) {
				//~ for (auto &atom3 : atom2) {
					//~ if (&atom3 != &atom1) {
						//~ for (auto &atom4 : atom2) { // impropers
							//~ if (&atom4 != &atom3 && &atom4 != &atom1) {
								//~ if (visited_impropers.count(make_tuple(&atom3, &atom1, &atom2, &atom4)) 
									//~ || visited_impropers.count(make_tuple(&atom1, &atom3, &atom2, &atom4)) 
									//~ || visited_impropers.count(make_tuple(&atom1, &atom4, &atom2, &atom3)) 
									//~ || visited_impropers.count(make_tuple(&atom4, &atom1, &atom2, &atom3)) 
									//~ || visited_impropers.count(make_tuple(&atom4, &atom3, &atom2, &atom1)) 
									//~ || visited_impropers.count(make_tuple(&atom3, &atom4, &atom2, &atom1))) continue;
								//~ visited_impropers.insert(make_tuple(&atom4, &atom3, &atom2, &atom1));
//~ 
								//~ this->improper.push_back(make_tuple(&atom3, &atom1, &atom2, &atom4));
								//~ dbgmsg("info.improper (improper) (third atom is central) = " 
									//~ << atom3.atom_number() << " " << atom1.atom_number() 
									//~ << " " << atom2.atom_number() << " " << atom4.atom_number());
							//~ }
						//~ }
					//~ }
				//~ }
			//~ }
		//~ }
		dbgmsg("BONDS ARE : " << endl << this->bond);
		dbgmsg("ANGLES ARE : " << endl << this->angle);
		dbgmsg("DIHEDRALS ARE : " << endl << this->dihedral);
		dbgmsg("IMPROPERS ARE : " << endl << this->improper);
		return *this;
	}
	//~ MoleculeInfo& MoleculeInfo::get_kb_force_info(const Molib::Molecule &receptor, 
					//~ const Molib::Molecule &ligand, const ForceField &ff,
					//~ const int &dist_cutoff) {
	MoleculeInfo& MoleculeInfo::get_kb_force_info(const Molib::Molecule &receptor, 
					const Molib::Molecule &ligand, const int &dist_cutoff) {
		dbgmsg("get_kb_force_info");
		// 1-2 and 1-3 bonded exclusions for knowledge-based forcefield
		Molib::AtomVec rec_atoms = receptor.get_atoms();
		Molib::AtomVec lig_atoms = ligand.get_atoms();
		Molib::AtomVec atoms;
		atoms.reserve(rec_atoms.size() + lig_atoms.size());
//~ #ifndef NDEBUG
		//~ for (auto &pa : rec_atoms) {
			//~ Molib::Atom &atom1 = *pa;
			//~ dbgmsg("receptor atom1 = " << atom1);
		//~ }
		//~ for (auto &pa : lig_atoms) {
			//~ Molib::Atom &atom1 = *pa;
			//~ dbgmsg("ligand atom1 = " << atom1);
		//~ }
//~ #endif
		dbgmsg("receptor atoms are : " << endl << rec_atoms);
		dbgmsg("ligand atoms are : " << endl << lig_atoms);

		atoms.insert(atoms.end(), rec_atoms.begin(), rec_atoms.end());
		atoms.insert(atoms.end(), lig_atoms.begin(), lig_atoms.end());
		//~ dbgmsg("get_kb_force_info2");
		//~ set<pair<Molib::Atom*, Molib::Atom*>> bonded_exclusions;
		set<Molib::AtomPair> bonded_exclusions;
		for (auto &pa : atoms) {
			Molib::Atom &atom1 = *pa;
			//~ dbgmsg("get_kb_force_info atom1 = " << atom1);
			for (auto &atom2 : atom1) {
			//~ for (auto &bond12 : atom1) {
				//~ auto &atom2 = bond12.second_atom();
				bonded_exclusions.insert({&atom1, &atom2});
				for (auto &atom3 : atom2) {
				//~ for (auto &bond23 : atom2) {
					//~ auto &atom3 = bond23.second_atom();
					bonded_exclusions.insert({&atom1, &atom3});
				}
			}
		}
//~ #ifndef NDEBUG
		//~ dbgmsg("get_kb_force_info3");
		//~ // calculate max distance to compute forces
		//~ const int sz = ((ff.kb_force_type.begin()->second).begin()->second).potential.size(); // number 
		//~ const double max_dist = ff.step * sz;
		//~ dbgmsg("max distance to compute knowledge-based ff " << max_dist 
				//~ << " currently selected distance cutoff is " << dist_cutoff);
//~ #endif
		//~ // knowledge-based forces between atoms except between bonded exclusions
		//~ Molib::MolGrid g(atoms);
		//~ for (auto &pa1 : atoms) {
			//~ for (auto &pa2 : g.get_neighbors(*pa1, dist_cutoff)) {
				//~ if (!bonded_exclusions.count({pa1, pa2})) {
					//~ this->kbforce.push_back({pa1, pa2});
				//~ }
			//~ }
		//~ }
		// knowledge-based forces between atoms except between bonded exclusions
		set<Molib::AtomPair> visited;
		Molib::MolGrid g(atoms);
		for (auto &pa1 : atoms) {
			for (auto &pa2 : g.get_neighbors(*pa1, dist_cutoff)) {
				if (!visited.count({pa1, pa2})) {
					visited.insert({pa2, pa1});
					if (!bonded_exclusions.count({pa1, pa2})) {
						this->kbforce.push_back({pa1, pa2});
					}
				}
			}
		}
		dbgmsg("KBFORCES ARE : " << endl << this->kbforce);
		return *this;
	}
	MoleculeInfo& MoleculeInfo::get_interaction_force_info(const Molib::Molecule &receptor, 
		const Molib::Molecule &ligand, const int &dist_cutoff) {

		// knowledge-based forces only between ligand and receptor atoms
		Molib::MolGrid g(receptor.get_atoms());
		for (auto &pa1 : ligand.get_atoms()) {
			for (auto &pa2 : g.get_neighbors(*pa1, dist_cutoff)) {
				this->kbforce.push_back({pa1, pa2});
			}
		}
		return *this;
	}
};
