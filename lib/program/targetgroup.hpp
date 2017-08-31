#ifndef TARGETGROUP_H
#define TARGETGROUP_H

#include "molib/it.hpp"
#include "program/target.hpp"

namespace Program {
        class CANDOCK_EXPORT TargetGroup {
                std::vector< std::unique_ptr<Target> > __targets;
        public:
                TargetGroup (const std::string &input_name);

                void dock_fragments(const FragmentLigands& ligands);
                void dock_ligands(const FragmentLigands& ligands);

                std::multiset<std::string>   determine_overlapping_seeds (const int max_seeds, const int number_of_occurances) const;
                std::set<std::string>        determine_non_overlapping_seeds (const TargetGroup &antitargets);
                void make_scaffolds(const TargetGroup& antitargets);
                void design_ligands(const TargetGroup& antitargets);
        };
}

#endif
