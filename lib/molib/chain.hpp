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

#ifndef CHAIN_H
#define CHAIN_H
#include "geom3d/geom3d.hpp"
#include "geom3d/matrix.hpp"
#include "it.hpp"
#include "atom.hpp"
#include "residue.hpp"

namespace Molib {
	class Residue;
	class Atom;
	class Model;
	
	class Chain : public template_map_container<Residue, Chain, Model, Residue::res_pair> {
		char __chain_id;
		Geom3D::Coordinate __crd; // geometric center
	public:
		Chain(char chain_id) : __chain_id(chain_id) {}
		Chain(const Chain &rhs) : __chain_id(rhs.__chain_id), __crd(rhs.__crd) { 
			for (auto &residue : rhs) { 
				dbgmsg("Copy constructor : chain");
				add(new Residue(residue)); 
			} 
		}
		void init_bio(Chain &chain_asym, const Geom3D::Matrix &bio_rota);
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false);
		
		Geom3D::Coordinate& crd() { return __crd; }
		Residue& add(Residue *r) { return this->aadd(Residue::res_pair(r->resi(), r->ins_code()), r, this); }
		Residue& residue(Residue::res_pair p) const { return this->element(p); }

		void set_crd(); // calculate geom center

		bool has_residue(Residue::res_pair p) { return this->has_element(p); }
		char chain_id() const { return __chain_id; }
		Atom::Vec get_atoms(const Residue::res_type &rest=Residue::res_type::notassigned) const;
		Chain& erase_properties() { for (auto &residue : *this) residue.erase_properties(); return *this; }
		
		friend ostream& operator<< (ostream& stream, const Chain& c);
	};

} // Molib
#endif
	
