#ifndef _PROBE_H
#define _PROBE_H

#include "sphere.h"
#include "const.h"


#ifdef CILE
class Probe;
//class Grid;
//class Molecule;
class EElement;
bool by_probe_dist(const Probe*, const Probe*);
#endif

class Probe : public Sphere {
 public:
 Probe(int NUM) : Sphere(NUM), next(NULL), size_move(0), color(0) {}
  Probe *next;
  int size_move;
  int color;
  Probe *move[5];
  void reset_visited();
  void output(int=0, bool=false);
  void invert(Probe*&);
  void output_probe();
  void print_probe(ostream &os=cout);
  void free();
#ifdef CILE
  EElement *neighb;
  float dist;
  Probe *previous;
  bool inQ;
  bool too_close;
//  void init_grid(Grid*, Molecule*, float);
  void pdb();
//  void delete_neighbor_list();
#endif
};

extern Probe *probe1, *probe2;


#endif // _PROBE_H
