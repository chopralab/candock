#ifndef SMILES_H
#define SMILES_H

#include <ostream>
#include <string>
#include <vector>

namespace candock {

namespace help {
struct edge {
    std::string atom_property1;
    std::string atom_property2;
    std::string bond_property;
    // string bond_stereo;
};
typedef std::vector<edge> smiles;
struct rename_rule {
    smiles pattern;
    std::vector<std::string> rule;
};
typedef std::vector<rename_rule> rename_rules;

std::ostream& operator<<(std::ostream& os, const smiles& edges);
std::ostream& operator<<(std::ostream& os, const rename_rule& rule);
}
}

#endif
