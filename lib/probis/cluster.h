#ifndef _CLUSTER_H
#define _CLUSTER_H

#include "const.h"

class Subgraph;
class Item;

class Cluster {
 public:
  void cluster(Subgraph*);
 private:
  int num_common(Item*, Item*);
};

extern Cluster *cluster;



#endif // _CLUSTER_H
