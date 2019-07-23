/* This is targetgroup.hpp and is part of CANDOCK
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

#ifndef TARGETGROUP_H
#define TARGETGROUP_H

#include "statchem/molib/it.hpp"
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
