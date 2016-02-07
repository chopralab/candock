#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include "nrset.hpp"
#include "hydrogens.hpp"
#include "geom3d/geom3d.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
#include "model.hpp"
#include "assembly.hpp"
#include "molecule.hpp"
#include "molecules.hpp"

using namespace std;

namespace Molib {

	ostream& operator<< (ostream& stream, const NRset& m) {
		for (auto &molecules : m) {
			stream << molecules;
		}
		return stream;
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


};
