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
