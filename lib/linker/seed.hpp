/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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
#include "helper/debug.hpp"
#include "molib/it.hpp"
#include "fragmenter/fragmenter.hpp"
#include "geom3d/coordinate.hpp"
#include "molib/internal.hpp"
#include "molib/molecule.hpp"
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
