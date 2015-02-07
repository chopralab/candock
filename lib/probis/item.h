#ifndef _ITEM_H
#define _ITEM_H

#include "const.h"
#include "vert.h"

class Molecule;

class Item {
 public:
// Item(int CLID, int SC, double BS, float RMS, float CS, Molecule *PM, gsl_matrix *ROTA, gsl_vector *T) : cluster_id(CLID), scons(SC), blosum_score(BS), calpha_rmsd(RMS), cluster_score(CS), pm(PM), U(ROTA), t(T) {}
  Item() { U = gsl_matrix_alloc(3, 3); t = gsl_vector_alloc(3); }

  ~Item() { gsl_matrix_free(U); gsl_vector_free(t);
            for (unsigned int i = 0; i < vert.size(); i++) delete vert[i]; 
            vert.clear();
            edge.clear();
            cons.clear();
          }

  vector< Vert* > vert;  // vert[vi]; vi ... vertex of this subgraph
  vector<pair<int,int > > edge; // edges in this subgraph (vi, vj)
  vector<pair<ChResi, int> > cons; // se uporablja v get_rota_trans_resi in output_rotate_pdb

  float calpha_rmsd, 
    surf_vector_angle, 
    cluster_score;

  double blosum_score;  // mora biti vecja natancnost (double), ker so lahko zelo majhne vrednosti

  float score_descriptor_ratio, score_probe;  // neuporabljeno!
  int N;        // iz kolikih klastrov je ta klaster
  gsl_matrix *U;
  gsl_vector *t;
  void copy_unique(Item*);
};




#endif // _ITEM_H

