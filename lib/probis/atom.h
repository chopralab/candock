#ifndef _ATOM_H
#define _ATOM_H

#include "const.h"
#include "sphere.h"
#include "colors.h"

class Coor;
class EElement;

class Atom : public Sphere {
//  class Model {
//    int id;
//    string pdb_id;
//    string chain_id;
//  public:
//  Model(int N, string P, string C) : num(N), pdb_id(P), chain_id(C) {}
//  } m;

 public:
 Atom(int NUM) : Sphere(NUM), next(NULL), neighb(NULL), acid_start(NULL), het(false), model(1) {}
// Atom(int NUM) : Sphere(NUM), m(1, "", ""), next(NULL), neighb(NULL), acid_start(NULL), het(false) {}
  ~Atom();
  Atom *next;
  EElement *neighb;
  Colors colors;
  char resn[4], tag[4], chain_id;
  int resi, aasize; 
  Atom *acid_start;
  bool het;
  int model;

  Atom  & operator = (const Atom &);
  bool overlap(Coor);
  void delete_neighbor_list();
  void print_atom(ostream &os=cout);
  void check_aasize();
#ifdef CILE
  void pdb();
#endif

};

#endif // _ATOM_H
