#ifndef SEED_H
#define SEED_H
#include "helper/debug.hpp"
#include "pdbreader/it.hpp"
#include "fragmenter/fragmenter.hpp"
#include "geom3d/coordinate.hpp"
#include "pdbreader/internal.hpp"
#include "pdbreader/molecule.hpp"
#include "graph/graph.hpp"
#include "segment.hpp"
#include <tuple>
#include <functional>

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
#endif
