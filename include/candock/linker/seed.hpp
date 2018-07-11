#ifndef SEED_H
#define SEED_H
#include <functional>
#include <tuple>
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/geometry/coordinate.hpp"
#include "candock/graph/graph.hpp"
#include "candock/helper/debug.hpp"
#include "candock/linker/segment.hpp"
#include "candock/molib/internal.hpp"
#include "candock/molib/it.hpp"
#include "candock/molib/molecule.hpp"

namespace candock {

namespace linker {
class State;
class Segment;

class Seed : public molib::template_vector_container<Seed*, Seed> {
    Segment& __seg;

   public:
    typedef graph::Graph<Seed> Graph;

    Seed(Segment& seg) : __seg(seg) {}
    Segment& get_segment() const { return __seg; }
    const std::string get_label() const {
        return __seg.get_label();
    }                                 // graph ostream operator
    int weight() const { return 0; }  // dummy for graph ostream operator
    friend std::ostream& operator<<(std::ostream& stream, const Seed& s);
    static Graph create_graph(const Segment::Graph& segment_graph);
};
};
}

#endif
