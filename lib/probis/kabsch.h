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

