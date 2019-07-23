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

#ifndef _PRODUCT_H
#define _PRODUCT_H


#include "const.h"

class Row;
class Item;
class Score;
class Descriptor;


/* napovemo primerjalne funkcije */
bool by_mnsp (Descriptor*, Descriptor*);
bool by_desc_num (Descriptor*, Descriptor*);




class Column {
 public:
  Column () : item(NULL), next(NULL) {}
  int num;
  Row *item;
  Column *next;
};

class Oglisca {
 public:
  Descriptor *d1, *d2, *d3, *d4;
  float a, b, c;   // stranice dolzina(a) > dolzina(b) > dolzina(c)
  float V; // volumen tetraedra
};


class Row {
 public:
 Row() : desc1(NULL), desc2(NULL), column(NULL), next(NULL) {}
  double score, combined_score;  // OBSOLETE
  int sgraph;
  int num;
  double surf_rmsd;
  int debug;
  Descriptor *desc1, *desc2;
  Column *column;
  Row *next;
};

class Product {
 public:
  Product() : row(NULL) {}
  ~Product();
  Row *row;
  int num_subgraph;
  int size;
  void tetrahedron(Score*, Descriptor*, Descriptor*, int);
//  void free();
  void count_subgraphs();
  void init_product(Score*);
 private:
  multimap<long int,Oglisca> thedron1;
  multimap<long int,Oglisca> thedron2;
  vector<pair<pair<Oglisca,Oglisca>, Item* > > filtrat;
  
};

extern Product *product;




#endif // _PRODUCT_H

