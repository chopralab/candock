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

#ifndef _CLIQUE_H
#define _CLIQUE_H


#include "const.h"
#include "bit.h"
#include "kabsch.h"

class Vert;
class Subgraph;


bool by_degree (Vert*, Vert*);


class Clique : public Kabsch, public Bit {
 public:
  Clique () { wipe_connect(); }
  ~Clique () { }
  void read_dimacs(Subgraph*, char [], int);
  void max_clique(Subgraph*);
  void MCQ(Item*, Item*);
  void score_descriptor_ratio(Subgraph*, Descriptor*, Descriptor*);
 private:
  void EXPAND(vector<Vert* > &, int);
  void COLOR_SORT(vector<Vert* > & );
  int CUT1(int, vector< Vert* > &);
  int CUT2(int, vector<Vert* > &, vector<Vert* > &);
  void COPY(vector<Vert* > &, vector<Vert* > & );
  void DELETE(vector<Vert* > & );
  void clear_connect(Item*);
  void wipe_connect();
  void assign_connect(Item*);
  Item *Q, *QMAX;
  int pk;

  clock_t start;
  bool end_condition;
};

extern Clique *clique;


#endif // _CLIQUE_H
