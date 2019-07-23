/* This is seed.hpp and is part of CANDOCK
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

#ifndef SEED_H
#define SEED_H

#include <functional>
#include <tuple>

#include "statchem/fragmenter/fragmenter.hpp"
#include "statchem/geometry/coordinate.hpp"
#include "statchem/graph/graph.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/molib/internal.hpp"
#include "statchem/molib/it.hpp"
#include "statchem/molib/molecule.hpp"

#include "candock/linker/segment.hpp"

namespace candock {
namespace linker {

using namespace statchem;

class State;
class Segment;

class Seed : public molib::template_vector_container<Seed*, Seed> {
    Segment& __seg;

   public:
    typedef statchem::graph::Graph<Seed> Graph;

    Seed(Segment& seg) : __seg(seg) {}
    Segment& get_segment() const { return __seg; }
    const std::string get_label() const {
        return __seg.get_label();
    }                                 // graph ostream operator
    int weight() const { return 0; }  // dummy for graph ostream operator
    friend std::ostream& operator<<(std::ostream& stream, const Seed& s);
    static Graph create_graph(const Segment::Graph& segment_graph);
};

}
}

#endif
