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

