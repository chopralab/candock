/* This is dockfragments.hpp and is part of CANDOCK
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

#ifndef DOCKFRAGMENTS_H
#define DOCKFRAGMENTS_H

#include "candock/program/programstep.hpp"

#include "candock/program/findcentroids.hpp"
#include "candock/program/fragmentligands.hpp"

#include "candock/docker/dock.hpp"
#include "statchem/molib/nrset.hpp"
#include "statchem/score/score.hpp"

namespace candock {

namespace Program {

class DockFragments : public ProgramStep {
    const FindCentroids& __found_centroids;
    const FragmentLigands& __fragmented_ligands;

    const score::Score& __score;
    const molib::Atom::Grid& __gridrec;

    const std::string& __name;
    std::string __top_seeds_location;

    molib::NRset __all_seeds;

    void __dock_fragment(int start, const docker::Gpoints& gpoints,
                         const docker::Gpoints& gpoints0);

   protected:
    virtual bool __can_read_from_files();
    virtual void __read_from_files();
    virtual void __continue_from_prev();

   public:
    DockFragments(const FindCentroids& found_centroids,
                  const FragmentLigands& fragmented_ligands,
                  const score::Score& score, const molib::Atom::Grid& gridrec,
                  const std::string& name);

    virtual ~DockFragments() {}

    std::vector<std::pair<double, std::string>> get_best_seeds() const;

    molib::NRset get_negative_seeds(const std::set<std::string>& seeds,
                                    const double max_value) const;
    molib::NRset get_top_seeds(const std::set<std::string>& seeds,
                               const double top_percent) const;
    molib::NRset get_seeds(const molib::Molecule& ligand,
                           const double top_percent) const;

    docker::Gpoints get_gridhcp();
};
}
}

#endif  // DOCKFRAGMENTSSTEP_H
