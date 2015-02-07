#ifndef _CLUSTERDATA_H
#define _CLUSTERDATA_H

#include "const.h"

class Molecule;

class ClusterData {
 public:
  ~ClusterData() { gsl_matrix_free(U); gsl_vector_free(t); }
 ClusterData() : cluster_id(0), scons(0), blosum_score(0), calpha_rmsd(0), surf_vector_angle(0), cluster_score(0), m(NULL), U(NULL), t(NULL) {}
 ClusterData(int CLID, int SC, double BS, float RMS, float SVA, float CS, Molecule *M, gsl_matrix *ROTA, gsl_vector *T) : cluster_id(CLID), scons(SC), blosum_score(BS), calpha_rmsd(RMS), surf_vector_angle(SVA), cluster_score(CS), m(M), U(ROTA), t(T) {}
  // ClusterData(int CLID, int SC, double BS, float RMS, float CS) : cluster_id(CLID), scons(SC), blosum_score(BS), calpha_rmsd(RMS), cluster_score(CS) {}
  int cluster_id, scons;
  double blosum_score;
  float calpha_rmsd;
  float surf_vector_angle;
  float cluster_score;
  Molecule *m;
  gsl_matrix *U;
  gsl_vector *t;
};

#endif // _CLUSTERDATA_H


