#ifndef SEED_H
#define SEED_H
#include "candock/helper/debug.hpp"
#include "candock/molib/it.hpp"
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/geometry/coordinate.hpp"
#include "candock/molib/internal.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/graph/graph.hpp"
#include "candock/linker/segment.hpp"
#include <tuple>
#include <functional>

namespace candock {

namespace Linker {
	class State;
	class Segment;
	
	class Seed : public template_vector_container<Seed*, Seed> {
		Segment &__seg;
	public:

		typedef Glib::Graph<Seed> Graph;

		Seed(Segment &seg) : __seg(seg) {}
		Segment& get_segment() const { return __seg; } 
		const string get_label() const { return __seg.get_label(); } // graph ostream operator
		int weight() const { return 0; } // dummy for graph ostream operator
		friend ostream& operator<< (ostream& stream, const Seed& s);
		static Graph create_graph(const Segment::Graph &segment_graph);
	
	};
};

}

#endif
