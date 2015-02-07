#ifndef INTERNAL_H
#define INTERNAL_H
#include "helper/error.hpp"
#include "helper/help.hpp"
#include "helper/debug.hpp"
#include "geom3d/geom3d.hpp"
#include "geom3d/coordinate.hpp"
#include <memory>
#include <queue>
#include <algorithm>


namespace Molib {
	class Atom;
	typedef map<const Atom*, Geom3D::Coordinate> AtomToCrd;
	struct Torsion {
		const Atom *a1, *a2, *a3, *a4;
		Torsion(const Atom *at1, const Atom *at2, const Atom *at3, const Atom *at4) : a1(at1), a2(at2), a3(at3), a4(at4) {}
		Torsion(const Torsion &other) { a1 = other.a1; a2 = other.a2; a3 = other.a3; a4 = other.a4; }
		Torsion& operator=(const Torsion &rv) { a1 = rv.a1; a2 = rv.a2; a3 = rv.a3; a4 = rv.a4; }
		bool operator==(const Torsion& right) const { return a1 == right.a1 && a2 == right.a2 && a3 == right.a3 && a4 == right.a4; }
		bool operator!=(const Torsion& right) const { return !(*this == right); }
	};
	extern const Torsion empty_torsion;
	class Internal {
		map<const Atom*, map<const Atom*, double>> __ic_bond;
		map<const Atom*, map<const Atom*, map<const Atom*, double>>> __ic_angle;
		map<const Atom*, map<const Atom*, map<const Atom*, map<const Atom*, double>>>> __ic_dihedral;
		Geom3D::Coordinate __set_crd(const Geom3D::Point&, const Geom3D::Point&, const Geom3D::Point&, const double, const double, const double) const;

	public:
		Internal() {}
		//~ Internal(const MolGraph &graph) { build(graph); }
		//~ void build(const MolGraph &);
		Internal(const AtomVec &atoms) { build(atoms); }
		void build(const AtomVec &atoms);
		AtomToCrd cartesian(const Atom&, const Atom&, const Atom&, const Geom3D::Coordinate&, const Geom3D::Coordinate&, const Geom3D::Coordinate&, const AtomSet&) const;
		void set_dihedral(const Atom &a1, const Atom &a2, const Atom &a3, const Atom &a4, const double angle) { __ic_dihedral[&a1][&a2][&a3][&a4] = angle; }
		double get_dihedral(const Atom &a1, const Atom &a2, const Atom &a3, const Atom &a4) const { return __ic_dihedral.at(&a1).at(&a2).at(&a3).at(&a4); }
		void set_dihedral(const Torsion &t, const double angle) { __ic_dihedral[t.a1][t.a2][t.a3][t.a4] = angle; }
		double get_dihedral(const Torsion &t) const { return __ic_dihedral.at(t.a1).at(t.a2).at(t.a3).at(t.a4); }
		//~ friend ostream& operator<< (ostream& stream, Graph<P>& g);
	};
};
#endif
