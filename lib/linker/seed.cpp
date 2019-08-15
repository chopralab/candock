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

#include "seed.hpp"
#include "segment.hpp"
#include "molib/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"

namespace Linker {
	ostream& operator<< (ostream& stream, const Seed& s) {
		return stream << "Seed = " << s.get_segment().get_seed_id() << endl;
	}

	Seed::Graph Seed::create_graph(const Segment::Graph &segment_graph) {
		dbgmsg("Create seed graph ...");
		vector<unique_ptr<Seed>> vertices;
		for (auto &s : segment_graph) { // add vertices to seed graph
			if (s.is_seed()) {
				vertices.push_back(unique_ptr<Seed>(new Seed(s)));
			}
		}				
		for (size_t i = 0; i < vertices.size(); ++i) { // ... edges ...
			for (size_t j = i + 1; j < vertices.size(); ++j) {
				const Segment &seg_i = vertices[i]->get_segment();
				const Segment &seg_j = vertices[j]->get_segment();
				dbgmsg("i = " << i << " j = " << j << " "
					<< seg_i << " " << seg_j << " is_seed_adjacent = "
					<< boolalpha << seg_i.is_seed_adjacent(seg_j));
				if (seg_i.is_seed_adjacent(seg_j)) {
					Seed &seed1 = *vertices[i];
					Seed &seed2 = *vertices[j];
					seed1.add(&seed2);
					seed2.add(&seed1);
				}
			}
		}
		return Seed::Graph(std::move(vertices), true, false);
	}


};
