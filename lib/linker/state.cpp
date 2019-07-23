/* This is state.cpp and is part of CANDOCK
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

#include "candock/linker/state.hpp"
#include "statchem/helper/benchmark.hpp"
#include "statchem/helper/help.hpp"
#include "candock/linker/segment.hpp"
#include "statchem/molib/molecule.hpp"
using namespace std;

namespace candock {
namespace linker {

State::Vec operator-(const State::Set& left, const State::Set& right) {
    State::Vec ret;
    for (auto& state : left)
        if (!right.count(state)) ret.push_back(state);
    return ret;
}

string State::pdb() const {
    stringstream ss;
    for (size_t i = 0; i < get_segment().get_atoms().size(); ++i) {
        molib::Atom& a = const_cast<molib::Atom&>(get_segment().get_atom(i));
        a.set_crd(get_crd(i));
        ss << a;
    }
    return ss.str();
}

bool State::clashes(const State& other, const double clash_coeff)
    const {  // clashes between this and other state
    const size_t other_size = other.get_crds().size();
    for (size_t i = 0; i < __crds.size(); ++i) {
        const geometry::Point& crd1 = get_crd(i);
        const molib::Atom& a1 = __segment.get_atom(i);
        ;
        if (__segment.is_common_atom(i)) {
            dbgmsg("common state atom = " << a1.atom_number()
                                          << " is not checked for clashes");
            continue;
        }
        const double vdw1 = a1.radius();
        for (size_t j = 0; j < other_size; ++j) {
            const geometry::Point& crd2 = other.get_crd(j);
            const molib::Atom& a2 = other.get_segment().get_atom(j);
            if (other.get_segment().is_common_atom(j)) {
                dbgmsg("common state atom = " << a2.atom_number()
                                              << " is not checked for clashes");
                continue;
            }
            const double vdw2 = a2.radius();
            if (crd1.distance_sq(crd2) < pow(clash_coeff * (vdw1 + vdw2), 2))
                return true;
        }
    }
    return false;
}

ostream& operator<<(ostream& stream, const State& s) {
    //~ stream << "State(address = " << &s <<", segment = " <<
    //s.__segment.get_seed_id() << ") "
    stream << "State(id = " << s.get_id()
           << ", segment = " << s.__segment.get_seed_id() << ") "
           << " energy = " << setprecision(4) << fixed << s.__energy
           << " atom_crd =  ";
    for (size_t i = 0; i < s.get_crds().size(); ++i) {
        const geometry::Point& crd = s.get_crd(i);
        const molib::Atom& a = s.__segment.get_atom(i);
        ;
        stream << a.atom_number() << " -> " << crd << " ";
    }
    return stream;
}

ostream& operator<<(ostream& stream, const State::Vec& sv) {
    for (auto& state : sv) stream << "MEMBER STATE : " << *state << endl;
    return stream;
}
};
}
