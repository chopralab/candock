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

#include "output.h"
#include "desc.h"
#include "atom.h"
#include "item.h"
#include "product.h"
#include "score.h"
#include "molecule.h"
#include "eelement.h"
#include "desc.h"
#include "probe.h"


void Output::srf_file(const string &srf_file, Molecule *mol, Descriptor *desc, 
	Probe *probe, bool _motif) {
	stringstream ss;
	ss << "MOL> " << mol->get_lcolor() << endl;
	for (Atom *atm = mol->atom; atm != NULL; atm = atm->next) {
		atm->print_atom(ss);
	}
	Descriptor *d = desc;
	while (d) {
		d->print_desc(ss);
		d = d->next;
	}
	for (Probe *p = probe; p != NULL; p = p->next) {
		if (_motif && mol->get_lcolor() != p->color) continue;
		p->print_probe(ss);
	}
	ofstream output_file(srf_file);
	output_file << ss.str();
	output_file.close();
}

void Output::out_score(Item *item) {
  cout << "START> " << "SCORE" << endl;
  char buffer[100];
  sprintf(&buffer[0], "%s", "                                                                                            ");

  sprintf(buffer, "SIZE = %d", (int) item->vert.size() );
  cout << buffer << endl;
  sprintf(buffer, "DESCRIPTOR RMSD = %.3f", item->calpha_rmsd);
  cout << buffer << endl;
  sprintf(buffer, "PROBE RMSD = %.3f", item->score_probe);
  cout << buffer << endl;
  sprintf(buffer, "DESCRIPTOR RATIO = %.3f", item->score_descriptor_ratio);
  cout << buffer << endl;
  sprintf(buffer, "SURFACE VECTORS ANGLE = %.3f", item->surf_vector_angle);
  cout << buffer << endl;
  sprintf(buffer, "BLOSUM SCORE = %.3e", item->blosum_score);
  cout << buffer << endl;
  sprintf(buffer, "CLUSTER SCORE = %.3f", item->cluster_score);
  cout << buffer << endl;

  cout << "END> " << "SCORE" << endl;
}

void Output::out_rota(gsl_matrix *U) {
  cout << "START> " << "ROTA" << endl;
  char buffer[100];
  sprintf(&buffer[0], "%s", "                                                                                            ");
//  if (item->U)
  for (int i = 0; i < 3; i++) {
    sprintf(buffer, "%12f%12f%12f", gsl_matrix_get(U, i, 0), gsl_matrix_get(U, i, 1), gsl_matrix_get(U, i, 2));
    cout << buffer << endl;
  }
//  
//  else
//    for (int i = 0; i < 3; i++) {
//      sprintf(buffer, "%12f%12f%12f", 0.0, 0.0, 0.0);
//      cout << buffer << endl;
//    }

  cout << "END> " << "ROTA" << endl;
}

void Output::out_trans(gsl_vector *t) {
  cout << "START> " << "TRANS" << endl;
  char buffer[100];
//  if (item->t) {
  sprintf(buffer, "%12f%12f%12f", gsl_vector_get(t, 0), gsl_vector_get(t, 1), gsl_vector_get(t, 2));
//  }
//  else
//    sprintf(buffer, "%12f%12f%12f", 0.0, 0.0, 0.0);
  cout << buffer << endl;
  cout << "END> " << "TRANS" << endl;

}

//void Output::out_rota(Item *item) {
//  cout << "START> " << "ROTA" << endl;
//  char buffer[100];
//  sprintf(&buffer[0], "%s", "                                                                                            ");
////  if (item->U)
//  for (int i = 0; i < 3; i++) {
//    sprintf(buffer, "%12f%12f%12f", gsl_matrix_get(item->U, i, 0), gsl_matrix_get(item->U, i, 1), gsl_matrix_get(item->U, i, 2));
//    cout << buffer << endl;
//  }
////  
////  else
////    for (int i = 0; i < 3; i++) {
////      sprintf(buffer, "%12f%12f%12f", 0.0, 0.0, 0.0);
////      cout << buffer << endl;
////    }
//
//  cout << "END> " << "ROTA" << endl;
//}
//
//void Output::out_trans(Item *item) {
//  cout << "START> " << "TRANS" << endl;
//  char buffer[100];
////  if (item->t) {
//  sprintf(buffer, "%12f%12f%12f", gsl_vector_get(item->t, 0), gsl_vector_get(item->t, 1), gsl_vector_get(item->t, 2));
////  }
////  else
////    sprintf(buffer, "%12f%12f%12f", 0.0, 0.0, 0.0);
//  cout << buffer << endl;
//  cout << "END> " << "TRANS" << endl;
//
//}

void Output::out_desc(Item *item) {
  /*
    Output the two aligned sets of descriptors
    Input: Item Q or QMAX
  */
  gsl_vector *vec = gsl_vector_alloc(3);  
  gsl_vector *y = gsl_vector_alloc(3);  
  int k = 1;
  cout << "START> " << "DESCRIPTORS" << " <START" << endl;
  for (unsigned int i = 0; i<item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

    gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
    gsl_vector_add(y, item->t);
//    output_row(y, 0.0, k++, item->vert[i]->row->desc1->atom->tag, "TYR", item->vert[i]->row->desc1->atom->helix);
    output_row(y, 0.0, k++, item->vert[i]->row->desc1->atom->tag, "TYR", item->vert[i]->row->desc1->atom->model);
    
    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);
    
//    output_row(vec, 1.0, k++, item->vert[i]->row->desc2->atom->tag, "TYR", item->vert[i]->row->desc2->atom->helix);
    output_row(vec, 1.0, k++, item->vert[i]->row->desc2->atom->tag, "TYR", item->vert[i]->row->desc2->atom->model);
  }
  cout << "END" << endl;
  cout << "END> " << "DESCRIPTORS" << " <END" << endl;
  gsl_vector_free(vec);
  gsl_vector_free(y);
}

void Output::out_desc_quiet(Item *item) {
  /*
    Output the two aligned sets of descriptors
    Input: Item Q or QMAX
  */
  char line[100];
  cout << "START> " << "DESC"<< endl;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
//    sprintf(line, "%3s%2d%4s%2d", item->vert[i]->row->desc1->atom->tag, item->vert[i]->row->desc1->atom->helix, item->vert[i]->row->desc2->atom->tag, item->vert[i]->row->desc2->atom->helix);
    sprintf(line, "%3s%2d%4s%2d", item->vert[i]->row->desc1->atom->tag, item->vert[i]->row->desc1->atom->model, item->vert[i]->row->desc2->atom->tag, item->vert[i]->row->desc2->atom->model);
    cout << line << endl;
  }
  cout << "END> " << "DESC" << endl;
}

void Output::out_desc_xquiet(Item *item) {
  /*
    Output the two aligned sets of descriptors. Only desc->num are printed!
    Input: Item Q or QMAX
  */
  char line[100];
  cout << "START> " << "DESC"<< endl;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    sprintf(line, "%5d%5d", item->vert[i]->row->desc1->num, item->vert[i]->row->desc2->num);
    cout << line << endl;
  }
  cout << "END> " << "DESC" << endl;
}

void Output::out_resi(Item *item) {
  gsl_vector *vec = gsl_vector_alloc(3);  
  gsl_vector *y = gsl_vector_alloc(3);  
  int k = 1;
  cout << "START> " << "RESIDUES" << " <START" << endl;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    for (Atom *atm = item->vert[i]->row->desc1->atom->acid_start; atm != NULL && atm->resi == item->vert[i]->row->desc1->atom->acid_start->resi; atm = atm->next) {
      gsl_vector_set(vec, 0, atm->crd.x);
      gsl_vector_set(vec, 1, atm->crd.y);
      gsl_vector_set(vec, 2, atm->crd.z);

      gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
      gsl_vector_add(y, item->t);
      output_row(y, 2.0, k++, atm->tag, atm->resn, atm->resi);
    }
  }
  cout << "END" << endl;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    for (Atom *atm = item->vert[i]->row->desc2->atom->acid_start; atm != NULL && atm->resi == item->vert[i]->row->desc2->atom->acid_start->resi; atm = atm->next) {
      gsl_vector_set(vec, 0, atm->crd.x);
      gsl_vector_set(vec, 1, atm->crd.y);
      gsl_vector_set(vec, 2, atm->crd.z);
      
      output_row(vec, 3.0, k++, atm->tag, atm->resn, atm->resi);
    }
  }
  cout << "END" << endl;
  cout << "END> " << "RESIDUES" << " <END" << endl;
  gsl_vector_free(vec);
  gsl_vector_free(y);
}

void Output::out_resi_quiet(Item *item, Score *score) {
  /*
    If no --blosum modifier is given, write out zeros in the blosum score field.
  */
  char line[100];
  cout << "START> " << "RESI (aa1 aa2 blosum)" << endl;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
//    if (_blosum)
    sprintf(line, "%3c%3c%5d%4c%3c%5d%5d", item->vert[i]->row->desc1->atom->chain_id, one_letter_code(item->vert[i]->row->desc1->atom->resn), item->vert[i]->row->desc1->atom->resi, 
            item->vert[i]->row->desc2->atom->chain_id, one_letter_code(item->vert[i]->row->desc2->atom->resn), item->vert[i]->row->desc2->atom->resi, 
            score->score_blosum(one_letter_code(item->vert[i]->row->desc1->atom->resn), one_letter_code(item->vert[i]->row->desc2->atom->resn)));
//    else
//      sprintf(line, "%3c%5d%4c%5d%5d", one_letter_code(item->vert[i]->row->desc1->atom->resn), item->vert[i]->row->desc1->atom->resi, one_letter_code(item->vert[i]->row->desc2->atom->resn), item->vert[i]->row->desc2->atom->resi, 0);
    cout << line << endl;
  }
  cout << "END> " << "RESI" << endl;
}

void Output::output_row(gsl_vector *vec, double color, int num, int resi) {
  char line[300];
  sprintf(line, "%-6s%5d%2s%-4s%3s%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f", "HETATM", num, " ", "C", "DES", " ", resi, gsl_vector_get(vec, 0), gsl_vector_get(vec, 1), gsl_vector_get(vec, 2), 1.0, color);
  cout << line << endl;
}

void Output::output_row(gsl_vector *vec, double color, int num) {
  char line[300];
  sprintf(line, "%-6s%5d%2s%-4s%3s%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f", "HETATM", num, " ", "C", "DES", " ", 1, gsl_vector_get(vec, 0), gsl_vector_get(vec, 1), gsl_vector_get(vec, 2), 1.0, color);
  cout << line << endl;
}

void Output::output_row(gsl_vector *vec, double color, int num, const char *tag, const char *resn, int resi) {
  char line[300];
  sprintf(line, "%-6s%5d%2s%-4s%3s%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f", "ATOM", num, " ", tag, resn, " ", resi, gsl_vector_get(vec, 0), gsl_vector_get(vec, 1), gsl_vector_get(vec, 2), 1.0, color);
  cout << line << endl;
}

//~ double Output::surf_rmsd_for_atoms(Atom *atom, Item *item, Molecule *mol) {
  //~ double rez = 0.0;
  //~ for (unsigned int i = 0; i < item->vert.size(); i++) {
    //~ if (mol == mol1) {
      //~ if (atom->crd == item->vert[i]->row->desc1->crd) 
        //~ rez += 1/pow(1.0, 2) * item->vert[i]->row->surf_rmsd;
      //~ else
        //~ rez += 1/pow(dist(atom->crd, item->vert[i]->row->desc1->crd), 2) * item->vert[i]->row->surf_rmsd;
    //~ }
    //~ else if (mol == mol2) {
      //~ if (atom->crd == item->vert[i]->row->desc2->crd) 
        //~ rez += 1/pow(1.0, 2) * item->vert[i]->row->surf_rmsd;
      //~ else
        //~ rez += 1/pow(dist(atom->crd, item->vert[i]->row->desc2->crd), 2) * item->vert[i]->row->surf_rmsd;
    //~ }
  //~ }
  //~ return rez;
//~ }

void Output::output_align_protein(Molecule *mol1, Molecule *mol2, Item *item) {
  /*
    Output the two aligned proteins; the first one is rotated to fit the second one.
    The matrix U and  vector t are used, from the subgraph->qmax[i], pointer to which is in item.
  */
  gsl_vector *vec = gsl_vector_alloc(3);  
  gsl_vector *y = gsl_vector_alloc(3);  
  int i;
  i = 0;
  for (Atom *atm = mol1->atom; atm != NULL; atm = atm->next) {
    gsl_vector_set(vec, 0, atm->crd.x);
    gsl_vector_set(vec, 1, atm->crd.y);
    gsl_vector_set(vec, 2, atm->crd.z);

    gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
    gsl_vector_add(y, item->t);
    //    output_row(y, surf_rmsd_for_atoms(atm, item, mol1), ++i, atm->tag, atm->resn, atm->resi);
    output_row(y, 4.0, ++i, atm->tag, atm->resn, atm->resi);
  }
  i = 0;
  for (Atom *atm = mol2->atom; atm != NULL; atm = atm->next) {
    gsl_vector_set(vec, 0, atm->crd.x);
    gsl_vector_set(vec, 1, atm->crd.y);
    gsl_vector_set(vec, 2, atm->crd.z);

    //    output_row(vec, surf_rmsd_for_atoms(atm, item, mol2), ++i, atm->tag, atm->resn, atm->resi);
    output_row(vec, 5.0, ++i, atm->tag, atm->resn, atm->resi);
  }
  cout << "END" << endl;
  gsl_vector_free(vec);
  gsl_vector_free(y);
}

void Output::output_align2(Item *item) {
  /*
    Output the two aligned sets of descriptors together with the surface probe atom centers.
    Input: Item Q or QMAX
  */
  gsl_vector *vec = gsl_vector_alloc(3);  
  gsl_vector *y = gsl_vector_alloc(3);  
  int k = 1;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

    gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
    gsl_vector_add(y, item->t);
    output_row(y, item->vert[i]->row->surf_rmsd, k++);

    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

    output_row(vec, item->vert[i]->row->surf_rmsd, k++);

    for (Element *tmp=item->vert[i]->row->desc1->neighb; tmp!=NULL; tmp=tmp->next) {
      gsl_vector_set(vec, 0, tmp->item->crd.x);
      gsl_vector_set(vec, 1, tmp->item->crd.y);
      gsl_vector_set(vec, 2, tmp->item->crd.z);
      
      gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
      gsl_vector_add(y, item->t);
      output_row(y, 0.2, k++);
      
    }
    for (Element *tmp=item->vert[i]->row->desc2->neighb; tmp!=NULL; tmp=tmp->next) {
      gsl_vector_set(vec, 0, tmp->item->crd.x);
      gsl_vector_set(vec, 1, tmp->item->crd.y);
      gsl_vector_set(vec, 2, tmp->item->crd.z);
      
      output_row(vec, 0.8, k++);
      
    }
  }
  cout << "END" << endl;
  gsl_vector_free(vec);
  gsl_vector_free(y);
}

void Output::output_align3(Item *item) {
  /*
    Output the two aligned sets of descriptors together with their sum vectors.
    Input: Item Q or QMAX
  */
  gsl_vector *vec = gsl_vector_alloc(3);  
  gsl_vector *y = gsl_vector_alloc(3);  
  int k = 1;
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);
    
    gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
    gsl_vector_add(y, item->t);
    output_row(y, 0.0, k++);
    
    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);
    
    output_row(vec, 1.0, k++);
    
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->sum.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->sum.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->sum.z);
    
    gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
    gsl_vector_add(y, item->t);
    output_row(y, 0.2, k++);
    
    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->sum.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->sum.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->sum.z);
    
    output_row(vec, 0.8, k++);
    
  }
  cout << "END" << endl;
  gsl_vector_free(vec);
  gsl_vector_free(y);
}

void Output::output_align4(Item *item) {
  /*
    Output the two aligned sets of descriptors together with the correspondent amino acids.
    Input: Item Q or QMAX
  */
  gsl_vector *vec = gsl_vector_alloc(3);  
  gsl_vector *y = gsl_vector_alloc(3);  
  int k = 1;

  cout << "MODEL        1"<<endl;

  for (unsigned int i = 0; i < item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

    gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
    gsl_vector_add(y, item->t);
    output_row(y, 0.0, k++, item->vert[i]->row->desc1->atom->resi);

    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

    output_row(vec, 1.0, k++, item->vert[i]->row->desc2->atom->resi);
  }

  cout << "TER" << endl;
  cout << "ENDMDL" << endl;
  cout << "MODEL        2"<<endl;

  for (unsigned int i = 0; i < item->vert.size(); i++) {
    for (Atom *atm = item->vert[i]->row->desc1->atom->acid_start; atm != NULL && atm->resi == item->vert[i]->row->desc1->atom->acid_start->resi; atm = atm->next) {
      gsl_vector_set(vec, 0, atm->crd.x);
      gsl_vector_set(vec, 1, atm->crd.y);
      gsl_vector_set(vec, 2, atm->crd.z);

      gsl_blas_dgemv(CblasNoTrans, 1, item->U, vec, 0, y);
      gsl_vector_add(y, item->t);
      output_row(y, 2.0, k++, atm->tag, atm->resn, atm->resi);
    }
  }


  cout << "TER" << endl;
  cout << "ENDMDL" << endl;
  cout << "MODEL        3"<<endl;

  for (unsigned int i = 0; i < item->vert.size(); i++) {
    for (Atom *atm = item->vert[i]->row->desc2->atom->acid_start; atm != NULL && atm->resi == item->vert[i]->row->desc2->atom->acid_start->resi; atm = atm->next) {
      gsl_vector_set(vec, 0, atm->crd.x);
      gsl_vector_set(vec, 1, atm->crd.y);
      gsl_vector_set(vec, 2, atm->crd.z);
      
      output_row(vec, 3.0, k++, atm->tag, atm->resn, atm->resi);
    }
  }

  cout << "TER" << endl;
  cout << "ENDMDL" << endl;
  cout << "END" << endl;

  gsl_vector_free(vec);
  gsl_vector_free(y);
}

