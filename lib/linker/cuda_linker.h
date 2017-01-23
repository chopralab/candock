#ifndef CUDA_LINKER_H
#define CUDA_LINKER_H

#include "linker.hpp"
#include "poses.hpp"
#include "geom3d/quaternion.hpp"
#include "score/score.hpp"
#include "pdbreader/nrset.hpp"
#include "pdbreader/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "graph/mcqd.hpp"
#include "modeler/modeler.hpp"
#include "geom3d/geom3d.hpp"
#include "cluster/greedy.hpp"
#include <queue>
#include "partial.hpp"

namespace Linker {

class cuda_linker {
public:
    
    void setup(const int segment_graph_size, const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter);
};

}
#endif
