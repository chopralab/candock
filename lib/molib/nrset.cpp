#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include "candock/molib/nrset.hpp"
#include "candock/geom3d/geom3d.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/residue.hpp"
#include "candock/molib/chain.hpp"
#include "candock/molib/model.hpp"
#include "candock/molib/assembly.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/molib/molecules.hpp"

using namespace std;

namespace Molib {

	ostream& operator<< (ostream& stream, const NRset& m) {
		for (auto &molecules : m) {
			stream << molecules;
		}
		return stream;
	}

	void NRset::jiggle(std::mt19937 rng) {
		std::uniform_real_distribution<> dist_thousandth( -1.0f / 1000.0f, 1.0f / 1000.0f);
		for (auto &patom : this->get_atoms()) {
			dbgmsg("before jiggle crd = " << patom->crd());
			Geom3D::Coordinate dcrd(dist_thousandth(rng), 
				dist_thousandth(rng),
				dist_thousandth(rng));
			patom->set_crd(patom->crd() + dcrd);
			dbgmsg("after jiggle crd = " << patom->crd());
		}
	}

	Molecule::Vec NRset::get_molecules(const Residue::res_type &rest) const {
		Molecule::Vec molecules;
		for (auto &mols : *this) {
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
