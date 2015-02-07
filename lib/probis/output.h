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
