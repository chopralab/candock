/* This is internal.hpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

#ifndef INTERNAL_H
#define INTERNAL_H
#include <algorithm>
#include <memory>
#include <queue>
#include "candock/geometry/coordinate.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/error.hpp"
#include "candock/molib/atom.hpp"

namespace candock {

namespace molib {

struct Torsion {
    const Atom *a1, *a2, *a3, *a4;
    Torsion(const Atom* at1, const Atom* at2, const Atom* at3, const Atom* at4)
        : a1(at1), a2(at2), a3(at3), a4(at4) {}
    Torsion(const Torsion& other) {
        a1 = other.a1;
        a2 = other.a2;
        a3 = other.a3;
        a4 = other.a4;
    }
    Torsion& operator=(const Torsion& rv) {
        a1 = rv.a1;
        a2 = rv.a2;
        a3 = rv.a3;
        a4 = rv.a4;
        return *this;
    }
    bool operator==(const Torsion& right) const {
        return a1 == right.a1 && a2 == right.a2 && a3 == right.a3 &&
               a4 == right.a4;
    }
    bool operator!=(const Torsion& right) const { return !(*this == right); }
};
extern const Torsion empty_torsion;
class Internal {
    std::map<const Atom*, std::map<const Atom*, double>> __ic_bond;
    std::map<const Atom*, std::map<const Atom*, std::map<const Atom*, double>>>
        __ic_angle;
    std::map<const Atom*,
             std::map<const Atom*,
                      std::map<const Atom*, std::map<const Atom*, double>>>>
        __ic_dihedral;
    geometry::Coordinate __set_crd(const geometry::Point&,
                                   const geometry::Point&,
                                   const geometry::Point&, const double,
                                   const double, const double) const;

    typedef std::map<const Atom*, geometry::Coordinate> AtomToCrd;

   public:
    Internal() {}
    Internal(const Atom::Vec& atoms) { build(atoms); }
    void build(const Atom::Vec& atoms);
    geometry::Point::Vec cartesian(const Atom& ini1, const Atom& ini2,
                                   const Atom& ini3,
                                   const geometry::Coordinate& crd1,
                                   const geometry::Coordinate& crd2,
                                   const geometry::Coordinate& crd3,
                                   const Atom::Vec& next_atoms) const;
    void set_dihedral(const Atom& a1, const Atom& a2, const Atom& a3,
                      const Atom& a4, const double angle) {
        __ic_dihedral[&a1][&a2][&a3][&a4] = angle;
    }
    double get_dihedral(const Atom& a1, const Atom& a2, const Atom& a3,
                        const Atom& a4) const {
        return __ic_dihedral.at(&a1).at(&a2).at(&a3).at(&a4);
    }
    void set_dihedral(const Torsion& t, const double angle) {
        __ic_dihedral[t.a1][t.a2][t.a3][t.a4] = angle;
    }
    double get_dihedral(const Torsion& t) const {
        return __ic_dihedral.at(t.a1).at(t.a2).at(t.a3).at(t.a4);
    }
};
};
}

#endif
