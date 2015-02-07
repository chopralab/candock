#ifndef _EELEMENT_H
#define _EELEMENT_H

#include "const.h"

class Sphere;
class Atom;
class Probe;


class Element {
 public:
 Element() : next(NULL), item(NULL) {}
  double dAB;
  Element *next;
  Sphere *item;
};

class EElement : public Element {
 public:
 EElement() : Element(), size(0) {}
  Atom *atm[15];
  Probe *probe[15];
  int size;
  void exchange(int, int);
};




#endif // _EELEMENT_H
