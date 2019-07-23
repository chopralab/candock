/* MIT License
*
* Copyright (c) 2017 Janez Konc
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

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
