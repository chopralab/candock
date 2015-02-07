#ifndef _DESC_H
#define _DESC_H

#include "sphere.h"
#include "schmitt.h"

class Column;
class Row;
class Coor;
class Element;
class Grid;
class Probe;

class Descriptor : public Sphere {
 public:
 Descriptor(int NUM) : Sphere(NUM), bb(false), s(new Schmitt()), column(NULL), neighb(NULL), next(NULL), current(NULL) {}
  ~Descriptor();
  int psurf, atom_num;
  bool bb;  // ali je backbone ?
  Schmitt *s;
  Column *column;
  Coor sum;

  Element *neighb;
  Atom *atom;
  Descriptor *next;
  Descriptor *ref;
  //  void output_descriptors(int);
  void output_descriptors();
  void init_grid_descriptors(Grid*, double, double);
  void init_grid_probe(Grid*, int, Probe*, double, double);
//  void set_schmitt_weights();
  void calculate_sum();
  void insert(Row*);
  void free();
//  void trim(string, string);
  Row* get();
  //  double compare(Descriptor*);
  void print_desc(ostream &os=cout);
  void output();
  void invert(Descriptor*&);
  void free_columns();
  void reset_visited();
  void mark_backbone();

//  void reset_current() {current = column; }
 private:
  Column *current;
  void free_neighb();
};

extern Descriptor *desc1, *desc2;


#endif // _DESC_H
