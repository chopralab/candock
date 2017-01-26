#include "cuda_linker.h"
#include "cuda.h"
#include "cuda_runtime.h"

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
#include "partial.hpp"
#include <queue>
#include <iostream>

namespace Linker{



void cuda_linker::setup(const int segment_graph_size, vector<unique_ptr<State>> &states, int iter, int num_states, int num_docked_seeds){
    //First load all data onto gpu

    //Partial *dev_start_conformation;
    //vector<unique_ptr<State>> *dev_states;
    //int dev_segment_graph_size, dev_inter;
    //cudaMalloc(&dev_states, sizeof(states) * num_states);
//    cudaMalloc(&dev_start_conformation, sizeof());
    //cudaMemcpy();
    
    
//    cudaFree(dev_states);

}
}
