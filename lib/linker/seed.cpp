#include "seed.hpp"
#include "segment.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"

namespace Molib {
	ostream& operator<< (ostream& stream, const Seed& s) {
		return stream << "Seed = " << s.get_segment().get_name() << endl;
	}

	Seed::Graph Seed::create_graph(const Segment::Graph &segment_graph, const Paths &paths) {
		dbgmsg("Create seed graph ...");
		vector<unique_ptr<Seed>> vertices;
		for (auto &s : segment_graph) // add vertices to seed graph
			if (s.is_seed()) 
				vertices.push_back(unique_ptr<Seed>(new Seed(s)));
				for (int i = 0; i < vertices.size(); ++i) { // ... edges ...
					for (int j = i + 1; j < vertices.size(); ++j) {
						const Segment &seg_i = vertices[i]->get_segment();
						const Segment &seg_j = vertices[j]->get_segment();
						dbgmsg("i = " << i << " j = " << j << " "
							<< seg_i << " " << seg_j << " path exists = "
							<< boolalpha << (paths.count({&seg_i, &seg_j})
							|| paths.count({&seg_j, &seg_i})));
						if (paths.count({&seg_i, &seg_j}) || paths.count({&seg_j, &seg_i})) {
#ifndef NDEBUG
						   Segment::Graph::Path path = paths.count({&seg_i, &seg_j}) ? 
							   paths.at({&seg_i, &seg_j}) : paths.at({&seg_j, &seg_i});
							dbgmsg("path = ");
							for (auto it = path.begin(); it != path.end(); ++it)
							dbgmsg(**it);
#endif
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
