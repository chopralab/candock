#ifndef CHAIN_H
#define CHAIN_H
#include "geom3d/geom3d.hpp"
#include "geom3d/matrix.hpp"
#include "it.hpp"
#include "atom.hpp"
#include "residue.hpp"

using namespace std;

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
	
