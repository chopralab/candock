/* This is poses.hpp and is part of CANDOCK
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

#ifndef POSES_H
#define POSES_H

#include "statchem/geometry/geometry.hpp"
#include "statchem/molib/grid.hpp"

#include "candock/linker/seed.hpp"
#include "candock/linker/segment.hpp"

namespace candock {

namespace linker {
class State;

class Poses {
    class AtomPoint {
        const geometry::Point& __crd;
        const molib::Atom& __atom;
        State& __state;

       public:
        AtomPoint(const geometry::Point& crd, const molib::Atom& atom,
                  State& state)
            : __crd(crd), __atom(atom), __state(state) {}
        const geometry::Point& crd() const { return __crd; }
        State& get_state() { return __state; }
        const molib::Atom& get_atom() const { return __atom; }
        void distance(double) const {}  // just dummy : needed by grid
        double radius() const { return __atom.radius(); }

        typedef std::vector<std::unique_ptr<AtomPoint>> UPVec;
        typedef std::vector<AtomPoint*> PVec;
        typedef statchem::molib::Grid<AtomPoint> Grid;
    };

    std::map<Segment::Id, AtomPoint::UPVec> __atompoints;
    std::map<Segment::Id, AtomPoint::Grid> __grid;

   public:
    Poses(const Seed::Graph& seed_graph);
    State::Set get_join_states(const State& state, Segment& segment2,
                               molib::Atom::Pair& jatoms,
                               const double max_linker_length,
                               const double lower_tol_seed_dist,
                               const double upper_tol_seed_dist);
};
}
}

#endif
