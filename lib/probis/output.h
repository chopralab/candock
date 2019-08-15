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

#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "const.h"

class Molecule;
class Item;
class Score;
class Atom;
class Descriptor;
class Probe;

class Output {
 public:
	void srf_file(const string&, Molecule*, Descriptor*, Probe*, bool);
  void output_align_protein(Molecule*, Molecule*, Item*);
  void output_align2(Item*);
  void output_align3(Item*);
  void output_align4(Item*);
  void out_desc(Item*);
  void out_resi(Item*);
  void out_desc_quiet(Item*);
  void out_desc_xquiet(Item*);
  void out_resi_quiet(Item*, Score*);
  void out_rota(gsl_matrix*);
  void out_trans(gsl_vector*);
  void out_score(Item*);
#ifdef VERB
  void output_row(gsl_vector*, double, int, int);
  void output_row(gsl_vector*, double, int);
#endif
 private:
#ifndef VERB
  void output_row(gsl_vector*, double, int, int);
  void output_row(gsl_vector*, double, int);
#endif
  void output_row(gsl_vector*, double, int, const char*, const char*, int);
  //~ double surf_rmsd_for_atoms(Atom*, Item*, Molecule*);
};


extern Output *out;


#endif // _OUTPUT_H
