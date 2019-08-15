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

#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
using namespace std;

namespace Molib {
	ostream& operator<< (ostream& stream, const Chain& c) {
		for (auto &residue : c) { 
			stream << residue;
		}
		stream << "TER" << endl;
		return stream;
	}

	void Chain::init_bio(Chain &chain_asym, const Geom3D::Matrix &bio_rota) {
		for (auto &residue : chain_asym) {
			Residue &last = add(new Residue(residue)); // pair and tuple !!
			last.init_bio(residue, bio_rota);
		}
	}
	
	void Chain::rotate(const Geom3D::Matrix &rota, const bool inverse) {
		for (auto &residue : *this) {
			residue.rotate(rota, inverse);
		}
	}

	void Chain::set_crd() {
		int sz=0; 
		for (auto &residue : *this) { 
			residue.set_crd(); 
			if (residue.rest() == Residue::protein || residue.rest() == Residue::nucleic) { 
				__crd = __crd + residue.crd(); 
				sz++; 
			}
		}
		if (sz !=0)
			__crd = __crd / sz;
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

};
