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
