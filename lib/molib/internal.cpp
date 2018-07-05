#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include "candock/geometry/geometry.hpp"
#include "candock/graph/graph.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/molib/internal.hpp"
using namespace std;

namespace candock {
namespace molib {
	const Torsion empty_torsion(nullptr, nullptr, nullptr, nullptr);

	void Internal::build(const Atom::Vec &atoms) {
		for (auto &pa1 : atoms) { // find all bonds lengths
			const Atom &a1 = *pa1;
			for (auto &a2 : a1) {
				__ic_bond[&a1][&a2] = a1.crd().distance(a2.crd());
				dbgmsg("bond length between atoms " << a1.atom_number() 
					<< " and " << a2.atom_number() << " is " << __ic_bond[&a1][&a2]);
				for (auto &a3 : a2) {  // find all angles
					if (&a1 != &a3) {
						__ic_angle[&a1][&a2][&a3] = geometry::angle(a1.crd(), a2.crd(), a3.crd());
						dbgmsg("angle between atoms " << a1.atom_number() << ", " 
							<< a2.atom_number() << " and " << a3.atom_number() << " is " 
							<< geometry::degrees(__ic_angle[&a1][&a2][&a3]));
						for (auto &a4 : a3) {  // find all torsional angles
							if (&a1 != &a4 && &a2 != &a4) {
								__ic_dihedral[&a1][&a2][&a3][&a4] = 
									geometry::dihedral(a1.crd(), a2.crd(), a3.crd(), a4.crd());
								dbgmsg("dihedral between atoms " << a1.atom_number() << ", " 
									<< a2.atom_number() << ", " << a3.atom_number() << " and " 
									<< a4.atom_number() << " is " 
									<< geometry::degrees(__ic_dihedral[&a1][&a2][&a3][&a4]));
							}
						}
					}
				}
			}
		}
	}

	geometry::Coordinate Internal::__set_crd(const geometry::Point &p1, const geometry::Point &p2, const geometry::Point &p3, const double bond, const double angle, const double dihedral) const {
		geometry::Vector3 v1 = p2 - p1;
		geometry::Vector3 v2 = p3 - p2;
		geometry::Vector3 n = geometry::Coordinate::cross(v1, v2).norm();
		geometry::Vector3 n2 = geometry::Coordinate::cross(v2, n).norm();
		geometry::Vector3 v2n = v2.norm();
		const double yy = bond * sin(angle);
		return p3 - v2n * (bond * cos(angle)) - n2 * (yy * cos (dihedral)) + n * (yy * sin(dihedral));
	}

	geometry::Point::Vec Internal::cartesian(const Atom &ini1, const Atom &ini2, const Atom &ini3, 
										const geometry::Coordinate &crd1, const geometry::Coordinate &crd2, 
										const geometry::Coordinate &crd3, const Atom::Vec &next_atoms) const {

		Atom::Set atoms(next_atoms.begin(), next_atoms.end());
		AtomToCrd new_crd {{&ini1, crd1}, {&ini2, crd2}, {&ini3, crd3}};

		queue<tuple<const Atom*, const Atom*, const Atom*>> triples;
		triples.push(make_tuple(&ini1, &ini2, &ini3));

		while (!triples.empty()) {
			const Atom *a1 = get<0>(triples.front());
			const Atom *a2 = get<1>(triples.front());
			const Atom *a3 = get<2>(triples.front());
			dbgmsg("cartesian a1 = " << *a1);
			dbgmsg("cartesian a2 = " << *a2);
			dbgmsg("cartesian a3 = " << *a3);
			triples.pop();
			if (__ic_dihedral.count(a1) && __ic_dihedral.at(a1).count(a2) 
				&& __ic_dihedral.at(a1).at(a2).count(a3)) {
				for (auto &kv : __ic_dihedral.at(a1).at(a2).at(a3)) {
					const Atom *a4 = kv.first;
					if (atoms.count(const_cast<Atom*>(a4))) {
						if (!new_crd.count(a4)) {
							new_crd.insert({a4, this->__set_crd(new_crd[a1], 
															new_crd[a2], 
															new_crd[a3], 
															__ic_bond.at(a3).at(a4), 
															__ic_angle.at(a2).at(a3).at(a4), 
															__ic_dihedral.at(a1).at(a2).at(a3).at(a4))
											});
							dbgmsg("coords of atom " << a4->atom_number() << " are " << a4->crd()
							<< " bond length between " << a3->atom_number() << " and " << a4->atom_number() 
							<< " is " << __ic_bond.at(a3).at(a4) << " angle between " << a2->atom_number() 
							<< " and " << a3->atom_number() << " and " << a4->atom_number() << " is " 
							<< __ic_angle.at(a2).at(a3).at(a4) << " dihedral between " << a1->atom_number() 
							<< " and " << a2->atom_number() << " and " << a3->atom_number() << " and " 
							<< a4->atom_number() << " is " << __ic_dihedral.at(a1).at(a2).at(a3).at(a4));

							triples.push(make_tuple(a2, a3, a4));

						}
					}
				}
			}
		}
		new_crd.erase(&ini1);

		// order coordinates according to next's segment atom order
		geometry::Point::Vec crds;
		for (size_t i = 0; i < next_atoms.size(); ++i) {
			crds.push_back(new_crd.at(next_atoms[i]));
		}
		return crds;
	}
};
}
