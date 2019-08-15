/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include "nrset.hpp"
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
