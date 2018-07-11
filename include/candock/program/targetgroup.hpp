#ifndef TARGETGROUP_H
#define TARGETGROUP_H

#include "candock/molib/it.hpp"
#include "candock/program/target.hpp"

namespace candock {

namespace Program {
class CANDOCK_EXPORT TargetGroup {
    std::vector<Target*> __targets;

   public:
    TargetGroup(const std::string& input_name);
    virtual ~TargetGroup();

    void dock_fragments(const FragmentLigands& ligands);
    void dock_ligands(const FragmentLigands& ligands);

    std::multiset<std::string> determine_overlapping_seeds(
        const int max_seeds, const int number_of_occurances) const;
    std::set<std::string> determine_non_overlapping_seeds(
        const TargetGroup& antitargets);
    void make_scaffolds(const TargetGroup& antitargets,
                        FragmentLigands& ligands);
    void design_ligands(const TargetGroup& antitargets,
                        FragmentLigands& ligands);
};
}
}

#endif
