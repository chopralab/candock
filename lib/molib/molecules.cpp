#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "atomtype.hpp"
#include "bondtype.hpp"
#include "bond.hpp"
#include "hydrogens.hpp"
#include "geom3d/geom3d.hpp"
#include "fragmenter/fragmenter.hpp"
#include "fragmenter/unique.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
#include "model.hpp"
#include "assembly.hpp"
#include "molecule.hpp"
#include "molecules.hpp"

using namespace std;

namespace Molib {

	ostream& operator<< (ostream& stream, const Molecules& m) {
		stream << setw(6) << left << "REMARK   4 NRPDB " << m.name() << endl;
		for (auto &molecule : m) {
			stream << molecule;
		}
		stream << setw(6) << left << "REMARK   4 END" << endl;
		return stream;
	}

	set<int> Molecules::get_idatm_types(set<int> previous) const {
		for (auto &molecule : *this)
		for (auto &pa : molecule.get_atoms()) {
			previous.insert(pa->idatm_type());
		}
		return previous;
	}

	void Molecules::rotate(const Geom3D::Matrix &rota, const bool inverse) {
		for (auto &molecule : *this) {
			molecule.rotate(rota, inverse);
		}
	}

	double Molecules::compute_max_radius() const { 
		Geom3D::Coordinate gc = compute_geometric_center();
		double max_radius = -HUGE_VAL;
		for (auto &molecule : *this) {
			Atom::Vec atoms = molecule.get_atoms();
			for (auto &patom : atoms)
				max_radius = max(max_radius, patom->crd().distance(gc));
		}
		return max_radius;
	}

	Geom3D::Coordinate Molecules::compute_geometric_center() const { 
		Geom3D::Coordinate center;
		for (auto &molecule : *this) {
			center = center + molecule.compute_geometric_center();
		}
		return center / this->size();
	}
	
	Molecules& Molecules::compute_idatm_type() { 
		auto &mols = *this;
		for (size_t i = 0; i < mols.size(); ++i) {
			auto &molecule = mols[i];
			try {
				AtomType::compute_idatm_type(molecule.get_atoms());
			} catch (exception &e) {
				log_error << "errmesg : deleting molecule " << molecule.name() 
					<< " (computing idatm types failed) due to " << e.what() << endl;
				mols.erase_shrink(i--);
			}
		}
		return *this; 
	}

	Molecules& Molecules::compute_hydrogen() { 
		auto &mols = *this;
		for (size_t i = 0; i < mols.size(); ++i) {
			auto &molecule = mols[i];
			try {
				for (auto &presidue : molecule.get_residues()) {
					Residue &residue = *presidue;
					if (!(help::standard_residues.count(residue.resn())
					   || help::ions.count(residue.resn()))) {
						Hydrogens::compute_hydrogen(residue);
					}
				}
			} catch (exception &e) {
				log_error << "errmesg : deleting molecule " << molecule.name() 
					<< " (computing hydrogens failed)" << endl;
				mols.erase_shrink(i--);
			}
		}
		return *this; 
	}

	Molecules& Molecules::compute_bond_order() { 
		auto &mols = *this;
		for (size_t i = 0; i < mols.size(); ++i) {
			auto &molecule = mols[i];
			try {
				for (auto &presidue : molecule.get_residues()) {
					Residue &residue = *presidue;
					if (!(help::standard_residues.count(residue.resn())
					   || help::ions.count(residue.resn()))) {
						BondOrder::compute_bond_order(residue.get_atoms(false));
					}
				}
			} catch (exception &e) {
				log_error << "errmesg : deleting molecule " << molecule.name() 
					<< " (bond order assignment failed) because " << e.what() << endl;
				mols.erase_shrink(i--);
			}
		}
		return *this; 
	}

	Molecules& Molecules::compute_bond_gaff_type() { 
		for (auto &presidue : this->get_residues()) {
			Residue &residue = *presidue;
			if (!(help::standard_residues.count(residue.resn())
			   || help::ions.count(residue.resn()))) {
				BondOrder::compute_bond_gaff_type(residue.get_atoms());
			}
		}
		return *this; 
	}

	Molecules& Molecules::refine_idatm_type() { 
		for (auto &presidue : this->get_residues()) {
			Residue &residue = *presidue;
			if (!(help::standard_residues.count(residue.resn())
			   || help::ions.count(residue.resn()))) {
				AtomType::refine_idatm_type(residue.get_atoms());
			}
		}
		return *this; 
	}

	Molecules& Molecules::erase_hydrogen() { 
		for (auto &presidue : this->get_residues()) {
			Residue &residue = *presidue;
			// Remove hydrogens added to protonated ions as well
			if (!(help::standard_residues.count(residue.resn())
			   || help::ions.count(residue.resn()))) {
				Hydrogens::erase_hydrogen(residue);
			}
		}
		return *this;
	}


	Molecules& Molecules::compute_ring_type() { 
		for (auto &presidue : this->get_residues()) {
			Residue &residue = *presidue;
			if (!(help::standard_residues.count(residue.resn())
			   || help::ions.count(residue.resn()))) {
				AtomType::compute_ring_type(residue.get_atoms());
			}
		}
		return *this;
	}

	Molecules& Molecules::compute_gaff_type() { 
		for (auto &presidue : this->get_residues()) {
			Residue &residue = *presidue;
			if (!(help::standard_residues.count(residue.resn())
			   || help::metals.count(residue.resn()))) {
				AtomType::compute_gaff_type(residue.get_atoms());
			}
		}
		return *this;
	}
        Molecules& Molecules::compute_rotatable_bonds() { 
                for (auto &presidue : this->get_residues()) {
                        Residue &residue = *presidue;
                        if (!(help::standard_residues.count(residue.resn())
                           || help::cofactor_residues.count(residue.resn())
                           || help::ions.count(residue.resn()))) {
                                BondOrder::compute_rotatable_bonds(residue.get_atoms());
                        }
                }
                return *this;
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

	Atom::Vec Molecules::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		Atom::Vec atoms;
		for (auto &molecule : *this) {
			auto ret = molecule.get_atoms(chain_ids, rest, model_number);
			atoms.insert(atoms.end(), ret.begin(), ret.end());
		}
		return atoms;
	}

	Residue::Vec Molecules::get_residues() const {
		Residue::Vec residues;
		for (auto &molecule : *this)
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			residues.push_back(&residue);
		}
		return residues;
	}

	Geom3D::Point::Vec Molecules::get_crds(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const { 

		Geom3D::Point::Vec crds;
		for (auto &patom : this->get_atoms(chain_ids, rest, model_number)) 
			crds.push_back(patom->crd());
		return crds;
	}

	Molecules& Molecules::compute_overlapping_rigid_segments(const string &seeds_file) { 
		Unique u(seeds_file); 
		for (auto &molecule : *this) 
			molecule.compute_overlapping_rigid_segments(u); 
		return *this; 
	}

        void create_mols_from_seeds(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
                for (auto &molecule : mols)
                for (auto &assembly : molecule)
                for (auto &model : assembly) {
                        for (auto &fragment : model.get_rigid()) { // iterate over seeds
                                if (fragment.is_seed()) {
                                        dbgmsg("considering to add " << fragment.get_seed_id());
                                        if (!added.count(fragment.get_seed_id())) { // take seeds that haven't been docked already
                                                dbgmsg("added " << fragment.get_seed_id());
                                                added.insert(fragment.get_seed_id());
                                                // add to new molecules
                                                Molib::Molecule &seed = seeds.add(new Molib::Molecule(std::to_string(fragment.get_seed_id())));
                                                Molib::Assembly &a = seed.add(new Molib::Assembly(0));
                                                Molib::Model &mod = a.add(new Molib::Model(1));
                                                Molib::Chain &c = mod.add(new Molib::Chain('X'));
                                                Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
                                                for (const Molib::Atom *atom : fragment.get_all()) {
                                                        r.add(new Molib::Atom(*atom));
                                                }
                                                seed.regenerate_bonds(molecule);
                                        }
                                }
                        }
                }
        }

};
