#ifndef _KABSCH_H
#define _KABSCH_H


#include "const.h"
#include "geo.h"

class Item;
class Descriptor;


class Kabsch {
 public:
  Kabsch () { 
    X = gsl_matrix_alloc(MAX_VERTICES, 3); Y = gsl_matrix_alloc(MAX_VERTICES, 3);
  }
  ~Kabsch () { gsl_matrix_free(X); gsl_matrix_free(Y); }
  double superimpose_calpha(Item*);
  double superimpose_calpha_single_desc(Item*);
  double superimpose_single_desc(Item*);
  double superimpose_only_calpha(Item*);
  double superimpose_calpha_cbeta(Item*);
  double superimpose_desc(Item*);
  double superimpose_desc2(Item*);
  Coor rotate_vector(Coor, gsl_matrix*, gsl_vector*);
  Coor rotate_vector_inv(Coor, gsl_matrix*, gsl_vector*);
  bool compare_matrices(gsl_matrix*, gsl_matrix*, gsl_vector*, gsl_vector*);
  void invert(gsl_matrix*, gsl_vector*);
  void output(gsl_matrix*, gsl_vector*);
  void set_rota_trans(gsl_matrix*, gsl_vector*, char**);
//  void set_rota_trans(gsl_matrix*, gsl_vector*, vector<string>);
  
 protected:
  Descriptor* align_descriptors(Descriptor*, Descriptor*, Item*);
  int kabsch(unsigned int, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_vector*, double*);
  gsl_matrix *X, *Y;

};



#endif // _KABSCH_H

