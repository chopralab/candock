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


