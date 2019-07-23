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

#ifndef _SCORE_H
#define _SCORE_H

#include "const.h"
#include "kabsch.h"

class Item;
class Molecule;
class Subgraph;
class Descriptor;
class Probe;

class Score : public Kabsch {
 public:
//  Score () { }
//  ~Score () { }
  void get_blosum();
  void get_internal_blosum();
  int score_blosum(unsigned char, unsigned char);
  void bit_score(vector <Item* > &);
  void E_score(Molecule*, Molecule*, vector <Item* > &);
  void cluster_score(const vector <Item* > &item);
  void score_descriptor_ratio(Subgraph*, Descriptor*, Descriptor*);
  void score_surface_vector(Subgraph*, Probe*, Probe*);
  void score_probe(Subgraph*, Probe*);
 private:
  short int blosum_matrix[256][256];


};

extern Score *score;


#endif // _SCORE_H
