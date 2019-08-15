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

#include "score.h"
#include "item.h"
#include "molecule.h"
#include "subgraph.h"
#include "desc.h"
#include "probe.h"
#include "atom.h"
#include "product.h"
#include "eelement.h"

void Score::get_internal_blosum() {
  /* 
     Preberi substitucijsko matriko, glede na BLOSUM parameters v parameters.inp.
     Mozne vrednosti za BLOSUM so: blosum62 ali blosum80
  */

  string *bl;

  if (strcmp(BLOSUM, "blosum80") == 0) {
    bl = blosum80_3;
  }
  else if (strcmp(BLOSUM, "blosum62") == 0) {
    bl = blosum62;
  }
  else {
    throw Err("Error (SCORE) : Wrong BLOSUM parameter, possible values are blosum80 or blosum62", 7);
  }
    

  unsigned char A[25];

  /*read the one letter codes of aacids */
  sscanf(bl[0].c_str(), "%*c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c", 
         &A[0], &A[1],&A[2],&A[3],&A[4],&A[5],&A[6],&A[7],&A[8],&A[9],&A[10],&A[11],&A[12],&A[13],
         &A[14],&A[15],&A[16],&A[17],&A[18],&A[19],&A[20],&A[21],&A[22],&A[23]);

  /*read the matrix */

  for (int i = 1; i < 25; i++) {
    unsigned char idx = static_cast<unsigned char> (bl[i].at(0));
    sscanf(bl[i].c_str(), "%*c %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd", 
           &blosum_matrix[idx][A[0]],  &blosum_matrix[idx][A[1]],  &blosum_matrix[idx][A[2]],
           &blosum_matrix[idx][A[3]],  &blosum_matrix[idx][A[4]],  &blosum_matrix[idx][A[5]],
           &blosum_matrix[idx][A[6]],  &blosum_matrix[idx][A[7]],  &blosum_matrix[idx][A[8]],
           &blosum_matrix[idx][A[9]],  &blosum_matrix[idx][A[10]], &blosum_matrix[idx][A[11]],
           &blosum_matrix[idx][A[12]], &blosum_matrix[idx][A[13]], &blosum_matrix[idx][A[14]],
           &blosum_matrix[idx][A[15]], &blosum_matrix[idx][A[16]], &blosum_matrix[idx][A[17]],
           &blosum_matrix[idx][A[18]], &blosum_matrix[idx][A[19]], &blosum_matrix[idx][A[20]],
           &blosum_matrix[idx][A[21]], &blosum_matrix[idx][A[22]], &blosum_matrix[idx][A[23]]);
//    cout << buffer[0] << " " << blosum_matrix[buffer[0]][A[0]] << " " << A[0] << " " << A[1] << " " << A[23] << endl;
  }
}

void Score::get_blosum() {
  /* 
     Read the substitution matrix.
  */

  ifstream blosum(BLOSUM);
  unsigned char buffer[256];
  unsigned char A[25];
  if (!blosum.is_open()) { 
    throw Err("Error (SCORE) : Cannot open blosum matrix file.", 8); 
  }
  /*read the comments first */
  while (!blosum.eof()) {
    blosum.getline ((char*)buffer,200);
    if (buffer[0] != '#') break;
  }
  /*read the one letter codes of aacids */
  sscanf((char*)buffer, "%*c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c", 
         &A[0], &A[1],&A[2],&A[3],&A[4],&A[5],&A[6],&A[7],&A[8],&A[9],&A[10],&A[11],&A[12],&A[13],
         &A[14],&A[15],&A[16],&A[17],&A[18],&A[19],&A[20],&A[21],&A[22],&A[23]);
  /*read the matrix */
  while (!blosum.eof()) {
    blosum.getline ((char*)buffer,200);
    if (!isalpha(buffer[0])) break;
    sscanf((char*)buffer, "%*c %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd", 
           &blosum_matrix[buffer[0]][A[0]], &blosum_matrix[buffer[0]][A[1]], &blosum_matrix[buffer[0]][A[2]],
           &blosum_matrix[buffer[0]][A[3]], &blosum_matrix[buffer[0]][A[4]], &blosum_matrix[buffer[0]][A[5]],
           &blosum_matrix[buffer[0]][A[6]], &blosum_matrix[buffer[0]][A[7]], &blosum_matrix[buffer[0]][A[8]],
           &blosum_matrix[buffer[0]][A[9]], &blosum_matrix[buffer[0]][A[10]], &blosum_matrix[buffer[0]][A[11]],
           &blosum_matrix[buffer[0]][A[12]], &blosum_matrix[buffer[0]][A[13]], &blosum_matrix[buffer[0]][A[14]],
           &blosum_matrix[buffer[0]][A[15]], &blosum_matrix[buffer[0]][A[16]], &blosum_matrix[buffer[0]][A[17]],
           &blosum_matrix[buffer[0]][A[18]], &blosum_matrix[buffer[0]][A[19]], &blosum_matrix[buffer[0]][A[20]],
           &blosum_matrix[buffer[0]][A[21]], &blosum_matrix[buffer[0]][A[22]], &blosum_matrix[buffer[0]][A[23]]);
//    cout << buffer[0] << " " << blosum_matrix[buffer[0]][A[0]] << " " << A[0] << " " << A[1] << " " << A[23] << endl;
  }
  blosum.close();
}

int Score::score_blosum(unsigned char amino1, unsigned char amino2) {
  /*
    Calculate blosum score for each matching aminoacid. Input are two aminoacids, one-letter codes,
    the output is blosum score (bit score). 
    NOTE: The output is not stored in any array. It's calculated on the fly.
  */
  return blosum_matrix[amino1][amino2];
}

void Score::bit_score(vector <Item* > &item) {
  /*
    Calculate the bit score, BS = (0.343 * sum_of_blosum(patch) - log(0.17)) / (log(2) * size(patch)),
    for a sequence of matched residues in clus. Numbers are valid for blosum80 matrix.
    Each maching pair of acids is counted only once in the score.
    NOTE: does not calculate bit score andymore, but only Score/size of a patch
  */
//  Item *item = subgraph->clus;
  float scbl;
  map<int,float> unique_patch;
  map<int,float>::iterator it;
  
  for (unsigned int i = 0; i < item.size(); i++) {
//    if (item[i].vert) {   // clus are not all initialized
      scbl = 0.0;
      unique_patch.clear();
      for (unsigned int j = 0; j < item[i]->vert.size(); j++) {
        char a1 = one_letter_code(item[i]->vert[j]->row->desc1->atom->resn);
        char a2 = one_letter_code(item[i]->vert[j]->row->desc2->atom->resn);
        unique_patch[item[i]->vert[j]->row->desc1->atom->resi] = score_blosum(a1, a2);
//        cout << "BS[" << one_letter_code(item[i].vert[j].row->desc1->atom->resn) << "," << one_letter_code(item[i].vert[j].row->desc2->atom->resn) << "] = " << score_blosum(one_letter_code(item[i].vert[j].row->desc1->atom->resn), one_letter_code(item[i].vert[j].row->desc2->atom->resn)) << endl;
      }
      for ( it=unique_patch.begin() ; it != unique_patch.end(); it++ ) {
//        cout << (*it).first << " => " << (*it).second << endl;
        scbl+= (*it).second;
      }
      item[i]->blosum_score = scbl / unique_patch.size();
//    }
  }
}

void Score::cluster_score(const vector <Item* > &item) {
  /* 
     Izracunamo cluster_score za vsak klaster (ali kliko), po enacbi :
     cluster_score = log(size * log( 1 + 1 / score_blosum ) * 1 / score_rmsd)
     size = stevilo verteksov v klastru
     Mora priti Score::E_score() 
   */
  for (size_t i = 0; i < item.size(); i++) {

    double bl = (item[i]->blosum_score < XS) ? XS : item[i]->blosum_score;  // da preprecimo nan
    float crms = (item[i]->calpha_rmsd < 1e-10) ? 1e-10 : item[i]->calpha_rmsd;  // da preprecimo nan

    item[i]->cluster_score = log ( item[i]->vert.size() * log( 1 + 1/bl ) * 1/crms );
  }
  
}

void Score::E_score(Molecule *m1, Molecule *m2, vector <Item* > &item) {
  /*
    Calculate the E score, E = K*m*n*exp(-lambda*S), where m and n are the lengths of the two compared sequences
    and S is the sum of patches' residue scores; K = 0.17 and lambda = 0.343 and  are valid for blosum80 matrix.
    Each matching pair of acids is counted only once in the score.
  */

  float S;
  map<ChResi,float> unique_patch;
  map<ChResi,float>::iterator it;

#ifdef VERB
  cout << "Score::E_score m = " << m1->sequence.size() << " n = " << m2->sequence.size() << endl;
#endif

  for (unsigned int i = 0; i < item.size(); i++) {

    S = 0.0;
    unique_patch.clear();

    for (unsigned int j = 0; j < item[i]->vert.size(); j++) {

      char a1 = one_letter_code(item[i]->vert[j]->row->desc1->atom->resn);
      char a2 = one_letter_code(item[i]->vert[j]->row->desc2->atom->resn);

      unique_patch[ChResi(item[i]->vert[j]->row->desc1->atom->chain_id, item[i]->vert[j]->row->desc1->atom->resi)] = score_blosum(a1, a2);
//      cout << " a1 = " << a1 << " a2 = " << a2 << endl;


    }


    for ( it=unique_patch.begin() ; it != unique_patch.end(); it++ ) {

      S += it->second;

    }

    item[i]->blosum_score = K * m1->sequence.size() * m2->sequence.size() * exp( (-1) * LAMBDA * S );

#ifdef VERB
    cout << "Score::E_score E-value = " << item[i]->blosum_score << " m1->seq = " << m1->sequence.size() << " m2->seq = " << m2->sequence.size() << endl;
#endif


  }
}


void Score::score_descriptor_ratio(Subgraph *subgraph, Descriptor *desc1, Descriptor *desc2) {
  /*
    Calculate and assign the ratio between matching and non-matching descriptors to each
    maximum clique. The matching and the non-matching descriptors N Angstroms from each 
    clique descriptor are counted. 
    OUTDATED NOTE: We compromise descriptor->psurf which is used to index different max cliques.
    NOTE: Now we use visited instead of psurf.
    WARNING: Has to be called before init_grid_probe and after init_grid_descriptors
  */
  Descriptor *d1, *d2;
  double num_matching, num_non_matching;
#ifdef VERB
  cout << "VERB> Score::score_descriptor_ratio subgraph->qmax.size() = " << subgraph->qmax.size() << endl;
#endif
  for (unsigned int i = 0; i < subgraph->qmax.size(); i++) {
#ifdef VERB
    cout << "VERB> Score::score_descriptor_ratio subgraph->qmax[" << i << "]->vert.size() = " << subgraph->qmax[i]->vert.size() << endl; 
#endif
    num_matching = 2 * subgraph->qmax[i]->vert.size();
    num_non_matching = 0;
    desc1->reset_visited();
    desc2->reset_visited();
    /* set 1 to all the descriptors visited in a qmax */
    for (unsigned int j = 0; j < subgraph->qmax[i]->vert.size(); j++) {
#ifdef VERB
      cout << "VERB> Score::score_descriptor_ratio j = " << j << endl;
#endif
      d1 = subgraph->qmax[i]->vert[j]->row->desc1;
      d2 = subgraph->qmax[i]->vert[j]->row->desc2;

      d1->visited = true;
      d2->visited = true;
    }
    /* calculate num_non_matching which are < CUTOFF_SECOND/2 away from clique descriptors */
    for (unsigned int j = 0; j < subgraph->qmax[i]->vert.size(); j++) {
      d1 = subgraph->qmax[i]->vert[j]->row->desc1;
      d2 = subgraph->qmax[i]->vert[j]->row->desc2;
      
      for (Element *tmp = d1->neighb; tmp!=NULL; tmp=tmp->next) {
        if ( ((Descriptor*)tmp->item)->visited == false 
             && dist(((Descriptor*)tmp->item)->crd, d1->crd) < CUTOFF_SECOND / 2) {
          num_non_matching++;
          ((Descriptor*)tmp->item)->visited = true;
        }
      }
      for (Element *tmp = d2->neighb; tmp!=NULL; tmp=tmp->next) {
        if ( ((Descriptor*)tmp->item)->visited == false
             && dist(((Descriptor*)tmp->item)->crd, d2->crd) < CUTOFF_SECOND / 2) {
          num_non_matching++;
          ((Descriptor*)tmp->item)->visited = true;
        }
      }
    }
    //    cout << "score=" << num_matching << endl;
    subgraph->qmax[i]->score_descriptor_ratio = (double) ( num_matching / ( num_matching + num_non_matching ));
  }
}

void Score::score_surface_vector(Subgraph *subgraph, Probe *probe1, Probe *probe2) {
  /*
    Calculate vector between average probe center and average descriptor center for each maximum clique found. 
    This function is intended to assess whether protein surfaces in the clique regions aren't superimposed 
    concave to convex.
   */
  Coor vector1, vector2;
  Descriptor *d1, *d2, *average_desc1, *average_desc2;
  Probe *average_probe1, *average_probe2;
  average_desc1 = new Descriptor(0);
  average_desc2 = new Descriptor(0);
  average_probe1 = new Probe(0);
  average_probe2 = new Probe(0);

  for (unsigned int i = 0; i < subgraph->qmax.size(); i++) {
    set_zero(average_probe1->crd);
    set_zero(average_probe2->crd);
    set_zero(average_desc1->crd);
    set_zero(average_desc2->crd);
    probe1->reset_visited();
    probe2->reset_visited();
    average_probe1->num = 0;
    average_probe2->num = 0;
    average_desc1->num = 0;
    average_desc2->num = 0;

    for (unsigned int j = 0; j < subgraph->qmax[i]->vert.size(); j++) {
      d1 = subgraph->qmax[i]->vert[j]->row->desc1;
      d2 = subgraph->qmax[i]->vert[j]->row->desc2;
      average_desc1->crd = average_desc1->crd + d1->crd;
      average_desc2->crd = average_desc2->crd + d2->crd;
      average_desc1->num++;
      average_desc2->num++;
      for (Element *tmp=d1->neighb; tmp!=NULL; tmp=tmp->next)
        if (!tmp->item->visited) {
          average_probe1->crd = average_probe1->crd + tmp->item->crd;
          average_probe1->num++;
          tmp->item->visited = true;
        }
      for (Element *tmp=d2->neighb; tmp!=NULL; tmp=tmp->next)
        if (!tmp->item->visited) {
          average_probe2->crd = average_probe2->crd + tmp->item->crd;
          average_probe2->num++;
          tmp->item->visited = true;
#ifdef VERB
          cout << "PRVIC: " << tmp->item->num << endl;
#endif
        }
#ifdef VERB
        else
          cout << "DVOJNIK: " << tmp->item->num << endl;
#endif
    }
    average_probe1->crd = average_probe1->crd * ((double) 1/average_probe1->num);  // scalar product
    average_probe2->crd = average_probe2->crd * ((double) 1/average_probe2->num);
    average_desc1->crd = average_desc1->crd * ((double) 1/average_desc1->num);
    average_desc2->crd = average_desc2->crd * ((double) 1/average_desc2->num);

    average_desc1->crd = rotate_vector(average_desc1->crd, subgraph->qmax[i]->U, subgraph->qmax[i]->t);
    average_probe1->crd = rotate_vector(average_probe1->crd, subgraph->qmax[i]->U, subgraph->qmax[i]->t);
    vector1 = average_probe1->crd - average_desc1->crd;
    vector2 = average_probe2->crd - average_desc2->crd;
    subgraph->qmax[i]->surf_vector_angle = angle(vector1, vector2, vector1 % vector2);
#ifdef VERB
    cout << "START> SURFACE VECTOR " << i << endl;
    Output *o = new Output();
    gsl_vector *v = gsl_vector_alloc(3);
    gsl_vector_set(v, 0, average_probe1->crd.x);
    gsl_vector_set(v, 1, average_probe1->crd.y);
    gsl_vector_set(v, 2, average_probe1->crd.z);
    o->output_row(v, 6.0, 1);
    gsl_vector_set(v, 0, average_desc1->crd.x);
    gsl_vector_set(v, 1, average_desc1->crd.y);
    gsl_vector_set(v, 2, average_desc1->crd.z);
    o->output_row(v, 6.0, 2);
    gsl_vector_set(v, 0, average_probe2->crd.x);
    gsl_vector_set(v, 1, average_probe2->crd.y);
    gsl_vector_set(v, 2, average_probe2->crd.z);
    o->output_row(v, 7.0, 3);
    gsl_vector_set(v, 0, average_desc2->crd.x);
    gsl_vector_set(v, 1, average_desc2->crd.y);
    gsl_vector_set(v, 2, average_desc2->crd.z);
    o->output_row(v, 7.0, 4);
    delete o;
    gsl_vector_free(v);
    cout << "END> SURFACE VECTOR" << endl;
#endif

#ifdef VERB
    cout << i << "-th SURFACE VECTOR ANGLES = " << subgraph->qmax[i]->score_surface_vector << endl;
    Coor t1, t2;
    t1.x = 1;t1.y=0;t1.z=0;
    t2.x = 1;t2.y=0;t2.z=0;
    cout << "ANGLE 0 = " << angle(t1, t2, t1 % t2) << endl;
    t1.x = 1;t1.y=1;t1.z=0;
    t2.x = 1;t2.y=0;t2.z=0;
    cout << "ANGLE 45 = " << angle(t1, t2, t1 % t2) << endl;
    t1.x = 0;t1.y=1;t1.z=0;
    t2.x = 1;t2.y=0;t2.z=0;
    cout << "ANGLE 90 = " << angle(t1, t2, t1 % t2) << endl;
    t1.x = -1;t1.y=1;t1.z=0;
    t2.x = 1;t2.y=0;t2.z=0;
    cout << "ANGLE 135 = " << angle(t1, t2, t1 % t2) << endl;
    t1.x = -1;t1.y=0;t1.z=0;
    t2.x = 1;t2.y=0;t2.z=0;
    cout << "ANGLE 180 = " << angle(t1, t2, t1 % t2) << endl;
    t1.x = -1;t1.y=-1;t1.z=0;
    t2.x = 1;t2.y=0;t2.z=0;
    cout << "ANGLE 225 = " << angle(t1, t2, t1 % t2) << endl;
#endif
    
  }
}

void Score::score_probe(Subgraph *subgraph, Probe *probe1) {
  /*
    Calculate the average RMSD between all probes in first and second protein belonging to the common regions. 
    NOTE: Only probe1 is needed, because we check if visited for it; but not for probe2.
  */
  Coor c1;
  int num;
  double d, min_dist;
  Descriptor *d1, *d2;

  for (unsigned int i = 0; i < subgraph->qmax.size(); i++) {
    probe1->reset_visited();
    //    probe2->reset_visited();
    num = 0;
    subgraph->qmax[i]->score_probe = 0.0;
    for (unsigned int j = 0; j < subgraph->qmax[i]->vert.size(); j++) {
      d1 = subgraph->qmax[i]->vert[j]->row->desc1;
      d2 = subgraph->qmax[i]->vert[j]->row->desc2;
      for (Element *tmp=d1->neighb; tmp!=NULL; tmp=tmp->next)
        if (!tmp->item->visited) {
          tmp->item->visited = true;
          c1 = rotate_vector(tmp->item->crd, subgraph->qmax[i]->U, subgraph->qmax[i]->t);
          min_dist = 100000.0;
          for (Element *tmp2=d2->neighb; tmp2!=NULL; tmp2=tmp2->next)
            if ((d = dist(c1, tmp2->item->crd)) < min_dist)
              min_dist = d;
#ifdef VERB
          cout << "score_probe: min_dist = " << min_dist << endl;
#endif
          if (min_dist != 100000) {
            subgraph->qmax[i]->score_probe += min_dist;
            num++;
          }
        }
    }
    if (num != 0)
      subgraph->qmax[i]->score_probe = subgraph->qmax[i]->score_probe / num;
    else 
      subgraph->qmax[i]->score_probe = 10.0;
  }
}
