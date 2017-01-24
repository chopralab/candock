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
