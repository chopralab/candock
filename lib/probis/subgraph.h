#ifndef _SUBGRAPH_H
#define _SUBGRAPH_H

#include "const.h"

class Output;
class Product;
class Score;
class Molecule;
class Item;
class Vert;
class Descriptor;

/* vsaka primerjalna (callback) funkcija mora biti napovedana */
bool by_hssp_length (const pair< TwoChResi, int>&, const pair< TwoChResi, int>&);
bool by_size (Item*, Item*);
bool by_cluster_score (Item*, Item*);
bool by_capacity (Item*, Item*);
bool by_size_blosum (Item*, Item*);
bool by_blosum (Item*, Item*);
bool by_size_surf_vect (Item*, Item*);
bool by_vertex (Vert*, Vert*);



/*
  INPUT: Graph G = (V,E) must given : V ... ARRAY[1..N] of vertices 
                                      E ... ARRAY of edges for each vertex (there are as many arrays as there are vertices)
*/

class Subgraph {
 public:
  Subgraph(Product*);
  Subgraph();
  ~Subgraph();
  vector<Item* > item, qmax, clus;
  void sort_qmax(bool (*)(Item*, Item*));
  void sort_clus(bool (*)(Item*, Item*));
  void sort_item(bool (*)(Item*, Item*));
  void output_qmax(Output*, Score* );
  void output_clus(Output*, Score* );
  void delete_bad_scoring();
  void extend_alignments(Score*, Molecule*, Molecule*, Product*, Descriptor*&, Descriptor*&, int);
  void extend_all(Score*, Molecule*, Molecule*, Product*, Descriptor*&, Descriptor*&);
  void extend_first_some(Score*, Molecule*, Molecule*, Product*, Descriptor*&, Descriptor*&);
};

extern Subgraph *subgraph;


#endif // _SUBGRAPH_H
