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
