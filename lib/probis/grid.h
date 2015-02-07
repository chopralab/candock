#ifndef _GRID_H
#define _GRID_H

#include "geo.h"

class Sphere;
class Probe;
class Element;

class Grid {
 public:
  Grid() : outsideGrid(false) {}
  //~ Grid (pair<Coor, Coor>);
  void init(pair<Coor, Coor>);
  ~Grid ();
  void deallocate_content();
  void make_grid(Sphere*);
  void volume_slice_desc(Sphere*, double, double);
  void volume_slice_atom(Sphere*, double, double);
  void volume_slice_surf(Probe*, double, double);
  void volume_slice_lig(Sphere*, double, double);
#ifdef CILE
  void volume_slice_probe(Probe*, double, double);
#endif
 private:

  typedef Element* pElement;

  pElement ***cell;
  Coor minCrd;
  Coor maxCrd;
  int cellX, cellY, cellZ;  // stevilo celic v smeri X, Y in Z (glej tudi SPACING)
  bool outsideGrid;

};

extern Grid *grid1, *grid2;

#endif // _GRID_H
