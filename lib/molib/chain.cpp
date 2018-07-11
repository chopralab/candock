#include "candock/molib/chain.hpp"
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "candock/molib/atom.hpp"
#include "candock/molib/residue.hpp"
using namespace std;

namespace candock {
namespace molib {
ostream& operator<<(ostream& stream, const Chain& c) {
    for (auto& residue : c) {
        stream << residue;
    }
    stream << "TER" << endl;
    return stream;
}

void Chain::init_bio(Chain& chain_asym, const geometry::Matrix& bio_rota) {
    for (auto& residue : chain_asym) {
        Residue& last = add(new Residue(residue));  // pair and tuple !!
        last.init_bio(residue, bio_rota);
    }
}

void Chain::rotate(const geometry::Matrix& rota, const bool inverse) {
    for (auto& residue : *this) {
        residue.rotate(rota, inverse);
    }
}

void Chain::set_crd() {
    int sz = 0;
    for (auto& residue : *this) {
        residue.set_crd();
        if (residue.rest() == Residue::protein ||
            residue.rest() == Residue::nucleic) {
            __crd = __crd + residue.crd();
            sz++;
        }
    }
    if (sz != 0) __crd = __crd / sz;
}

Atom::Vec Chain::get_atoms(const Residue::res_type& rest) const {
    Atom::Vec atoms;
    for (auto& residue : *this) {
        if (rest == Residue::res_type::notassigned || residue.rest() == rest) {
            auto ret = residue.get_atoms();
            atoms.insert(atoms.end(), ret.begin(), ret.end());
        }
    }
    return atoms;
}
};
}
