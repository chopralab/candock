/* This is linkfragments.hpp and is part of CANDOCK
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

#ifndef LINKFRAGMENTS_H
#define LINKFRAGMENTS_H

#include "statchem/modeler/forcefield.hpp"
#include "statchem/molib/molecule.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "candock/program/dockfragments.hpp"
#include "candock/program/programstep.hpp"
#include "statchem/score/score.hpp"

namespace candock {

namespace Program {

class LinkFragments : public ProgramStep {
    molib::Molecules __all_top_poses;
    const DockFragments& __seeds_database;

    const molib::Molecule& __receptor;
    const score::Score& __score;
    const OMMIface::ForceField& __ffield;
    const molib::Atom::Grid& __gridrec;

    int __ligand_cnt = 0;

    void __link_ligand(molib::Molecule& ligand);

   protected:
    virtual bool __can_read_from_files();
    virtual void __read_from_files();
    virtual void __continue_from_prev();

   public:
    LinkFragments(const molib::Molecule& receptor, const score::Score& score,
                  const OMMIface::ForceField& ffield,
                  const DockFragments& seeds_database,
                  const molib::Atom::Grid& gridrec)
        : __seeds_database(seeds_database),
          __receptor(receptor),
          __score(score),
          __ffield(ffield),
          __gridrec(gridrec){};

    virtual ~LinkFragments() {}
    void link_ligands(const molib::Molecules& ligands);
    const molib::Molecules& top_poses();
    void clear_top_poses() { __all_top_poses.clear(); };
};
}
}

#endif  // LINKFRAGMENTSSTEP_H
