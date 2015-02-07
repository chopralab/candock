#ifndef _SPHERE_H
#define _SPHERE_H

#include "geo.h"

class Sphere {
 public:
  Sphere (int NUM) : num(NUM), r(0.0), visited(false) {}
  Coor crd;
  int num;
  double r;
  bool visited;
 protected:
  void print_sphere(ostream &os=cout);
};

#endif // _SPHERE_H
