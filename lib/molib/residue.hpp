#ifndef RESIDUE_H
#define RESIDUE_H
#include "geom3d/geom3d.hpp"
#include "it.hpp"
#include "atom.hpp"

namespace Molib {
	class Atom;
	class Chain;
	
	class Residue : public template_map_container<Atom, Residue, Chain, int> {
	public:
		typedef enum {notassigned=1, protein=2, nucleic=4, ion=8, water=16, hetero=32} res_type;
		typedef tuple<char,std::string, int, char> res_tuple2;
		typedef pair<int, char> res_pair;
		typedef vector<Residue*> Vec;
		typedef set<Residue*> Set;
	private:
	std::string __resn;
		int __resi;
		char __ins_code;
		Residue::res_type __rest;
		Geom3D::Coordinate __crd; // geometric center
	public:
		Residue(string resn, int resi, char ins_code, res_type(rest)) 
			: __resn(resn), __resi(resi), __ins_code(ins_code), __rest(rest) {}
		Residue(const Residue &rhs) : __resn(rhs.__resn), __resi(rhs.__resi), 
			__ins_code(rhs.__ins_code), __rest(rhs.__rest), __crd(rhs.__crd) { 
			for (auto &atom : rhs) { 
				dbgmsg("Copy constructor : residue");
				add(new Atom(atom)); 
			} 
		}
		void init_bio(Residue &residue_asym, const Geom3D::Matrix &bio_rota);
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false);
		
		Geom3D::Coordinate& crd() { return __crd; }
		const Geom3D::Coordinate& crd() const { return __crd; }
		void set_crd() { for (auto &atom : *this) { __crd = __crd + atom.crd(); } 
			__crd = __crd / this->size(); } // calculate geom center
		Atom& add(Atom *a) { return this->aadd(a->atom_number(), a, this); }
	std::string resn() const { return __resn; }
		void set_resn(const std::string &resn) { __resn = resn; }
		void set_resi(int resi) { __resi = resi; }
		int resi() const { return __resi; }
		char ins_code() const { return __ins_code; }
		Residue::res_type rest() const { return __rest; }
		Atom& atom(int p) const { return this->element(p); }
		bool has_atom(int p) { return this->has_element(p); }
		Atom::Vec get_atoms(bool include_ions = true) const;
		void renumber_atoms(int new_start);
		Residue& erase_properties() { for (auto &atom : *this) atom.erase_properties(); return *this; }

		Residue& regenerate_bonds(const Residue&);

		friend ostream& operator<< (ostream& stream, const Residue& r);
	};


} // Molib
#endif
	
