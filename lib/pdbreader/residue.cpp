#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include "atomtype.hpp"
#include "bondtype.hpp"
#include "atom.hpp"
#include "residue.hpp"
using namespace std;

namespace Molib {

	ostream& operator<< (ostream& stream, const Residue& r) {
		for (auto &atom : r) { 
			stream << atom;
		} 
		return stream;
	}

	void Residue::init_bio(Residue &residue_asym, const Geom3D::Matrix &bio_rota) {
		for (Atom &atom : residue_asym) {
			Atom &last = add(new Atom(atom));
			last.crd().rotate_inline(bio_rota); // do the actual rotation
		}
	}
	void Residue::rotate(const Geom3D::Matrix &rota, const bool inverse) {
		if (inverse) for (auto &atom : *this) {	atom.crd().inverse_rotate_inline(rota); }
		else for (auto &atom : *this) {	atom.crd().rotate_inline(rota); }
	}

	Atom::Vec Residue::get_atoms() const {
		Atom::Vec atoms;
		for (auto &atom : *this) {
			atoms.push_back(&atom);
		}
		return atoms;
	}
};
