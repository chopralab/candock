/* This is state.hpp and is part of CANDOCK
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

#ifndef STATE_H
#define STATE_H

#include <functional>
#include <tuple>

#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/bond.hpp"

namespace candock {
namespace linker {

using namespace statchem;

class Segment;
class State {
   public:
    struct comp {
        bool operator()(State* const i, State* const j) const {
            return i->get_id() < j->get_id();
        }
    };

    typedef std::vector<State*> Vec;
    typedef std::vector<State> NVec;
    typedef std::set<State*, State::comp> Set;
    typedef std::vector<const State*> ConstVec;
    typedef std::pair<const State*, const State*> ConstPair;
    typedef int Id;

   private:
    const Segment& __segment;
    geometry::Point::Vec __crds;
    double __energy;
    Id __id;
#ifndef NDEBUG
    int __no;
#endif
   public:
    State(const Segment& segment, const geometry::Point::Vec crds,
          const double energy = 0)
        : __segment(segment), __crds(crds), __energy(energy) {}
    void set_energy(const double energy) { __energy = energy; }
    double get_energy() const { return __energy; }
    const Segment& get_segment() const { return __segment; }
    const geometry::Point::Vec& get_crds() const { return __crds; }
    geometry::Point::Vec& get_crds() { return __crds; }
    const geometry::Point& get_crd(const int i) const { return __crds[i]; }
    bool clashes(const State& other, const double clash_coeff)
        const;  // clashes between this and other state
    std::string pdb() const;
    void set_id(Id id) { __id = id; }
    Id get_id() const { return __id; }
#ifndef NDEBUG
    void set_no(int no) { __no = no; }
    int get_no() const { return __no; }
#endif
    friend std::ostream& operator<<(std::ostream& stream, const State& s);
    friend std::ostream& operator<<(std::ostream& stream, const Vec& sv);
};

State::Vec operator-(const State::Set& left, const State::Set& right);
};
}

#endif
