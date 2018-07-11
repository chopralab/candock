#include "candock/helper/smiles.hpp"

namespace candock {
namespace help {
std::ostream& operator<<(std::ostream& os, const smiles& edges) {
    for (auto& e : edges)
        os << '{' << e.atom_property1 << ',' << e.atom_property2 << ','
           << e.bond_property << '}';
    return os;
}

std::ostream& operator<<(std::ostream& os, const rename_rule& rule) {
    os << "PATTERN = " << rule.pattern << " RULE = {";
    for (auto& txt : rule.rule) os << txt << ',';
    os << '}';
    //~ os << endl;
    return os;
}
}
}
