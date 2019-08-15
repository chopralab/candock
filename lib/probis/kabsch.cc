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

#include "kabsch.h"
#include "desc.h"
#include "item.h"
#include "atom.h"
#include "product.h"


Descriptor* Kabsch::align_descriptors(Descriptor *desc1, Descriptor *desc2, Item *qmax) {
  /*
    Return the pointer to the two aligned sets of descriptors.
    Input: two descriptors sets and a maximum clique with the rotational and translational transformations
    Output: dtot
  */
  int i = 0;
  Descriptor *dtot = NULL, *tmp;
//  gsl_vector *vec = gsl_vector_alloc(3);  
//  gsl_vector *y = gsl_vector_alloc(3);  
  
  for (Descriptor *d1 = desc1; d1 != NULL; d1 = d1->next) {
//    gsl_vector_set(vec, 0, d1->crd.x);
//    gsl_vector_set(vec, 1, d1->crd.y);
//    gsl_vector_set(vec, 2, d1->crd.z);
//    
//    gsl_blas_dgemv(CblasNoTrans, 1, qmax->U, vec, 0, y);
//    gsl_vector_add(y, qmax->t);

    tmp = new Descriptor(i++);
    tmp->crd = rotate_vector(d1->crd, qmax->U, qmax->t);

//    tmp->crd.x = gsl_vector_get(y, 0);
//    tmp->crd.y = gsl_vector_get(y, 1);
//    tmp->crd.z = gsl_vector_get(y, 2);
    tmp->s->mnsp = d1->s->mnsp;
    //    cout << "mnsp = " << tmp->s->mnsp << "x = " << tmp->crd.x << " y = " << tmp->crd.y << " z = " << tmp->crd.y << endl;
    tmp->psurf = 0;
    tmp->ref = d1;
    tmp->next = dtot;
    dtot = tmp;
  }
  
  for (Descriptor *d2 = desc2; d2 != NULL; d2 = d2->next) {
    tmp = new Descriptor(i++);
    tmp->crd.x = d2->crd.x;
    tmp->crd.y = d2->crd.y;
    tmp->crd.z = d2->crd.z;
    tmp->s->mnsp = d2->s->mnsp;
    tmp->psurf = 1;
    tmp->ref = d2;
    tmp->next = dtot;
    dtot = tmp;
  }
  
//  gsl_vector_free(vec);
//  gsl_vector_free(y);
  return dtot;
}

/* gsl does not provide it */
static inline void gsl_vector_cross(
  const gsl_vector *a,
  const gsl_vector *b,
  gsl_vector *c
) {
  double a0=gsl_vector_get(a,0);
  double a1=gsl_vector_get(a,1);
  double a2=gsl_vector_get(a,2);
  double b0=gsl_vector_get(b,0);
  double b1=gsl_vector_get(b,1);
  double b2=gsl_vector_get(b,2);
  gsl_vector_set(c,0,a1*b2-b1*a2);
  gsl_vector_set(c,1,a2*b0-b2*a0);
  gsl_vector_set(c,2,a0*b1-b0*a1);
}

int Kabsch::kabsch(
  unsigned int size, /* the number of points */
  gsl_matrix *X,     /* the points to be moved */
  gsl_matrix *Y,     /* the points to move to */
  gsl_matrix *U,     /* the rotation matrix */
  gsl_vector *t,     /* the translation vector */
  double *s          /* the optimal scaling, if != 0 */
) {
  unsigned int i,j,k;
  int U_ok=1;
  double n=1.0/size;
  gsl_vector *cx=gsl_vector_alloc(3);     /* centroid of X */
  gsl_vector *cy=gsl_vector_alloc(3);     /* centroid of Y */
  gsl_matrix *R=gsl_matrix_alloc(3,3);    /* Kabsch's R */
  gsl_matrix *RTR=gsl_matrix_alloc(3,3);  /* R_trans * R (and Kabsch's bk) */
  gsl_eigen_symmv_workspace *espace=gsl_eigen_symmv_alloc(3);
  gsl_matrix *evec=gsl_matrix_alloc(3,3); /* eigenvectors (and Kabsch's ak) */
  gsl_vector *eval=gsl_vector_alloc(3);   /* vector of eigenvalues */
  /* compute centroid of X */
  gsl_vector_set_zero(cx);
  for(i=size;i>0;) {
    gsl_vector_const_view row=gsl_matrix_const_row(X,--i);
    gsl_vector_add(cx,&row.vector);
  } 
  gsl_vector_scale(cx,n);

  /* compute centroid of Y */
  gsl_vector_set_zero(cy);
  for(i=size;i>0;) {
    gsl_vector_const_view row=gsl_matrix_const_row(Y,--i);
    gsl_vector_add(cy,&row.vector);
  } 
  gsl_vector_scale(cy,n);


  /* move X to origin */
  for(i=size;i>0;) {
    gsl_vector_view row=gsl_matrix_row(X,--i);
    gsl_vector_sub(&row.vector,cx);
  }
  /* move Y to origin */
  for(i=size;i>0;) {
    gsl_vector_view row=gsl_matrix_row(Y,--i);
    gsl_vector_sub(&row.vector,cy);
  }

  if(size==1) {
    /* just one point, so U is trival */
    gsl_matrix_set_identity(U);
  }
  else {
    /* compute R */
    gsl_matrix_set_zero(R);
    for(k=size;k>0;) {
      --k;
      for(i=3;i>0;) {
        --i;
        for(j=3;j>0;) {
          --j;
          gsl_matrix_set(R,i,j,
            gsl_matrix_get(R,i,j)+
            gsl_matrix_get(Y,k,i)*gsl_matrix_get(X,k,j)
          );
        }
      }
    }

    /* compute RTR = R_trans * R */
    gsl_matrix_set_zero(RTR);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,R,0.0,RTR);

    /* compute orthonormal eigenvectors */
    gsl_eigen_symmv(RTR,eval,evec,espace);  /* RTR will be modified! */
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    if(gsl_vector_get(eval,1)>NORM_EPS) {
      /* compute ak's (as columns of evec) and bk's (as columns of RTR) */
      double norm_b0,norm_b1,norm_b2;
      gsl_vector_const_view a0=gsl_matrix_const_column(evec,0);
      gsl_vector_const_view a1=gsl_matrix_const_column(evec,1);
      gsl_vector_view a2=gsl_matrix_column(evec,2);
      gsl_vector_view b0=gsl_matrix_column(RTR,0);
      gsl_vector_view b1=gsl_matrix_column(RTR,1);
      gsl_vector_view b2=gsl_matrix_column(RTR,2);
      gsl_vector_cross(&a0.vector,&a1.vector,&a2.vector); /* a2 = a0 x a1 */
      gsl_blas_dgemv(CblasNoTrans,1.0,R,&a0.vector,0.0,&b0.vector);
      norm_b0=gsl_blas_dnrm2(&b0.vector);
      gsl_blas_dgemv(CblasNoTrans,1.0,R,&a1.vector,0.0,&b1.vector);
      norm_b1=gsl_blas_dnrm2(&b1.vector);
      if(norm_b0>NORM_EPS&&norm_b1>NORM_EPS) {
        gsl_vector_scale(&b0.vector,1.0/norm_b0);         /* b0 = ||R * a0|| */
        gsl_vector_scale(&b1.vector,1.0/norm_b1);         /* b1 = ||R * a1|| */
        gsl_vector_cross(&b0.vector,&b1.vector,&b2.vector);  /* b2 = b0 x b1 */

        norm_b2=gsl_blas_dnrm2(&b2.vector);
        if(norm_b2>NORM_EPS) {
          /* we reach this point only if all bk different from 0 */
          /* compute U = B * A_trans (use RTR as B and evec as A) */
          gsl_matrix_set_zero(U); /* to avoid nan */
          gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,RTR,evec,0.0,U);
        }
        else {
          U_ok=0;
          gsl_matrix_set_identity(U);
        }
      }
      else {
        U_ok=0;
        gsl_matrix_set_identity(U);
      }
    }
    else {
      U_ok=0;
      gsl_matrix_set_identity(U);
    }
  }

  /* compute t = cy - U * cx */
  gsl_vector_memcpy(t,cy);
  gsl_blas_dgemv(CblasNoTrans,-1.0,U,cx,1.0,t);


  if(s) {
    /* let us compute the optimal scaling as well */
    /* s = <Y,UX> / <UX,UX> */
    *s=1.0;
    if(U_ok&&size>1) {
      double dom=0.0;
      double nom=0.0;
      double dom_i,nom_i;
      gsl_vector *Uxi=gsl_vector_alloc(3);
      for(i=size;i>0;) {
        gsl_vector_const_view row_x=gsl_matrix_const_row(X,--i);
        gsl_vector_const_view row_y=gsl_matrix_const_row(Y,i);
        gsl_vector_set_zero(Uxi);
        gsl_blas_dgemv(CblasNoTrans,1.0,U,&row_x.vector,1.0,Uxi);
        gsl_blas_ddot(&row_y.vector,Uxi,&nom_i);
        nom+=nom_i;
        gsl_blas_ddot(Uxi,Uxi,&dom_i);
        dom+=dom_i;
      }
      *s=nom/dom;
      gsl_vector_free(Uxi);
    }
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  gsl_eigen_symmv_free(espace);
  gsl_matrix_free(RTR);
  gsl_matrix_free(R);
  gsl_vector_free(cy);
  gsl_vector_free(cx);

  return U_ok;
}

//double Kabsch::superimpose(Item *item) {
//  gsl_matrix *ROTA =  item->U; 
//  gsl_vector *trans = item->t;
//  gsl_vector *vec = gsl_vector_alloc(3);
//  for (unsigned int i = 0; i < item->vert.size(); i++) {
//    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
//    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
//    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);
//
//    gsl_matrix_set_row(X, i, vec);
//
//    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
//    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
//    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);
//
//#ifdef VERB
//    gsl_vector *temp = gsl_vector_alloc(3);
//
//    gsl_matrix_get_row(temp, X, i);
//    cout << "pred x: "<<  gsl_vector_get(temp, 0) << endl;
//
//    gsl_vector_free(temp);
//#endif
//
//    gsl_matrix_set_row(Y, i, vec);
//  }
//  gsl_vector_free(vec);
//  if (kabsch(item->vert.size(), X, Y, ROTA, trans, NULL) == 1) {
//
//#ifdef VERB
//    for (unsigned int i = 0; i < item->vert.size(); i++) {
//
//      gsl_vector *temp = gsl_vector_alloc(3);
//      
//      gsl_matrix_get_row(temp, X, i);
//      cout << "po x: "<<  gsl_vector_get(temp, 0) << endl;
//      
//      gsl_vector_free(temp);
//      
//    }
//#endif
//
//    return rmsd(item, ROTA, trans);
//
//  }
//  else 
//    return -1.0;
//
//}

double Kabsch::superimpose_desc(Item *item) {
  /*
    Prilegaj deskriptorje, vendar vrni rmsd med C-alpha atomi. Na ta nacin preprecis
    primere, ko sta prilegani funkcionalni skupini obrnjeni ena proti drugi.

    Upostevamo vec istih calpha atomov.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);

  /* najprej prilegamo deskriptorje */
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

    gsl_matrix_set_row(X, i, vec);

    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

    gsl_matrix_set_row(Y, i, vec);
  }


  /* ce prileganje uspe, potem superimposamo, po dobljenih matrikah, calpha atome, in izracunamo rmsd med calpha */
  if (kabsch(item->vert.size(), X, Y, ROTA, trans, NULL) == 1) {

    double rmsd = 0.0;

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      
      Coor c = rotate_vector(item->vert[i]->row->desc1->atom->acid_start->next->crd, ROTA, trans); 
      
      rmsd += dist_fast(c, item->vert[i]->row->desc2->atom->acid_start->next->crd);

    }

    gsl_vector_free(vec);

    /* vrnemo rmsd */
    return sqrt( rmsd  /  item->vert.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}

double Kabsch::superimpose_desc2(Item *item) {
  /*
    Prilegaj deskriptorje, vendar vrni rmsd med C-alpha atomi. Na ta nacin preprecis
    primere, ko sta prilegani funkcionalni skupini obrnjeni ena proti drugi.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);

  /* najprej prilegamo deskriptorje */
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

    gsl_matrix_set_row(X, i, vec);

    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

    gsl_matrix_set_row(Y, i, vec);
  }


  /* ce prileganje uspe, potem superimposamo, po dobljenih matrikah, calpha atome, in izracunamo rmsd med calpha */
  if (kabsch(item->vert.size(), X, Y, ROTA, trans, NULL) == 1) {

    set<pair<int, int> > unique_resi;
    unique_resi.clear();
  
    double rmsd = 0.0;

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
      int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;


      /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
      if (unique_resi.find( make_pair(resi1,resi2) ) == unique_resi.end() ) {

        unique_resi.insert(make_pair(resi1, resi2) );

        Coor c = rotate_vector(item->vert[i]->row->desc1->atom->acid_start->next->crd, ROTA, trans); 

        rmsd += dist_fast(c, item->vert[i]->row->desc2->atom->acid_start->next->crd);
      }

    }

    gsl_vector_free(vec);

    /* vrnemo rmsd */
    return sqrt( rmsd / unique_resi.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}


double Kabsch::superimpose_calpha(Item *item) {
  /*
    Funkcija prilega deskriptorje in calpha atome naenkrat, vrne pa rmsd
    med C-alpha atomi. Na ta nacin preprecis
    primere, ko sta prilegani funkcionalni skupini obrnjeni ena proti drugi.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);
  
  int st = 0;

  /* najprej dodamo deskriptorje v X in Y matriko */
  for (unsigned int i = 0; i < item->vert.size(); i++) {
    gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

    gsl_matrix_set_row(X, st, vec);

    gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
    gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
    gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

    gsl_matrix_set_row(Y, st, vec);
    
    st++;

  }

  /* 
     nato dodamo unique calpha v X in Y matriko, ce je na primer 2x ARG201--ARG345, dodamo samo enkrat calpha, 
     ce je ARG201--ARG345 in ARG201--ARG346, potem dodamo dvakrat calpha.

     NOVO APR/17/2010 : stevilka residueja ni vec edinstvena, ker imamo lahko vec chain-ov, sedaj razlikujemo residueje po acid_start
  */

  set<pair<Atom*, Atom*> > unique_resi;  // residueje razlikujemo po acid_start

  unique_resi.clear();

  for (unsigned int i = 0; i < item->vert.size(); i++) {

    Atom* acid_start1 = item->vert[i]->row->desc1->atom->acid_start;
    Atom* acid_start2 = item->vert[i]->row->desc2->atom->acid_start;

    /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
    if (unique_resi.find( make_pair(acid_start1,acid_start2) ) == unique_resi.end() ) {

      unique_resi.insert(make_pair(acid_start1, acid_start2) );
      gsl_vector_set(vec, 0, item->vert[i]->row->desc1->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc1->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc1->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(X, st, vec);
      
      gsl_vector_set(vec, 0, item->vert[i]->row->desc2->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc2->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc2->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(Y, st, vec);

      st++;

    }
  }

  /* ce prileganje uspe, alignamo calpha atome z izracunano matriko ROTA in trans in izracunamo RMSD med calpha */
  if (kabsch(st, X, Y, ROTA, trans, NULL) == 1) {

  
    double rmsd = 0.0;

    unique_resi.clear();

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      Atom* acid_start1 = item->vert[i]->row->desc1->atom->acid_start;
      Atom* acid_start2 = item->vert[i]->row->desc2->atom->acid_start;

      /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
      if (unique_resi.find( make_pair(acid_start1,acid_start2) ) == unique_resi.end() ) {


        unique_resi.insert(make_pair(acid_start1, acid_start2) );
      
        Coor c = rotate_vector(item->vert[i]->row->desc1->atom->acid_start->next->crd, ROTA, trans); 

        rmsd += dist_fast(c, item->vert[i]->row->desc2->atom->acid_start->next->crd);
      }
      
    }
    
    gsl_vector_free(vec);
    
    /* vrnemo rmsd */
    return sqrt( rmsd / unique_resi.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}

double Kabsch::superimpose_calpha_single_desc(Item *item) {
  /*
    Funkcija prilega deskriptorje (enega na vsak par prileganih ak, torej ce je
    vec, npr ARG13-ARG130 deskriptorjev vzame le enega izmed njih)
    in calpha atome naenkrat, vrne pa rmsd
    med C-alpha atomi. Na ta nacin preprecis
    primere, ko sta prilegani funkcionalni skupini obrnjeni ena proti drugi.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);
  
  int st = 0;

  set<pair<int, int> > unique_resi;
  unique_resi.clear();

  /* najprej dodamo deskriptorje v X in Y matriko */
  for (unsigned int i = 0; i < item->vert.size(); i++) {


    int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
    int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

    /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
    if (unique_resi.find( make_pair(resi1,resi2) ) == unique_resi.end() ) {

      unique_resi.insert(make_pair(resi1, resi2) );


      gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

      gsl_matrix_set_row(X, st, vec);

      gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

      gsl_matrix_set_row(Y, st, vec);
    
      st++;
    }

  }

  /* 
     nato dodamo unique calpha v X in Y matriko, ce je na primer 2x ARG201--ARG345, dodamo samo enkrat calpha, 
     ce je ARG201--ARG345 in ARG201--ARG346, potem dodamo dvakrat calpha.
  */

  unique_resi.clear();

  for (unsigned int i = 0; i < item->vert.size(); i++) {

    int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
    int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

    /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
    if (unique_resi.find( make_pair(resi1,resi2) ) == unique_resi.end() ) {

      unique_resi.insert(make_pair(resi1, resi2) );

      gsl_vector_set(vec, 0, item->vert[i]->row->desc1->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc1->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc1->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(X, st, vec);
      
      gsl_vector_set(vec, 0, item->vert[i]->row->desc2->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc2->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc2->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(Y, st, vec);

      st++;

    }
  }

  /* ce prileganje uspe, alignamo calpha atome z izracunano matriko ROTA in trans in izracunamo RMSD med calpha */
  if (kabsch(st, X, Y, ROTA, trans, NULL) == 1) {

  
    double rmsd = 0.0;

    unique_resi.clear();

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
      int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

      /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
      if (unique_resi.find( make_pair(resi1, resi2) ) == unique_resi.end() ) {

        unique_resi.insert( make_pair(resi1, resi2) );

        Coor c = rotate_vector(item->vert[i]->row->desc1->atom->acid_start->next->crd, ROTA, trans); 

        rmsd += dist_fast(c, item->vert[i]->row->desc2->atom->acid_start->next->crd);
      }
      
    }
    
    gsl_vector_free(vec);
    
    /* vrnemo rmsd */
    return sqrt( rmsd / unique_resi.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}


double Kabsch::superimpose_single_desc(Item *item) {
  /*
    Funkcija prilega deskriptorje (enega na vsak par prileganih ak, torej ce je
    vec, npr ARG13-ARG130 deskriptorjev vzame le enega izmed njih - na ta nacin 
    preprecimo preveliko utezenost ARG stranskih skupin.

    Vrne pa rmsd med deskriptorji.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);
  
  int st = 0;

  set<pair<int, int> > unique_resi;
  unique_resi.clear();

  /* najprej dodamo deskriptorje v X in Y matriko */
  for (unsigned int i = 0; i < item->vert.size(); i++) {


    int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
    int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

    /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
    if (unique_resi.find( make_pair(resi1,resi2) ) == unique_resi.end() ) {

      unique_resi.insert(make_pair(resi1, resi2) );


      gsl_vector_set(vec, 0, item->vert[i]->row->desc1->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc1->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc1->crd.z);

      gsl_matrix_set_row(X, st, vec);

      gsl_vector_set(vec, 0, item->vert[i]->row->desc2->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc2->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc2->crd.z);

      gsl_matrix_set_row(Y, st, vec);
    
      st++;
    }

  }
  /* ce prileganje uspe, alignamo deskriptorje z izracunano matriko ROTA in trans in izracunamo RMSD med deskriptorji */
  if (kabsch(st, X, Y, ROTA, trans, NULL) == 1) {

  
    double rmsd = 0.0;

    unique_resi.clear();

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
      int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

      /* upostevamo le en par deskriptorjev na vsak par prileganih AK */
      if (unique_resi.find( make_pair(resi1, resi2) ) == unique_resi.end() ) {

        unique_resi.insert( make_pair(resi1, resi2) );

        Coor c = rotate_vector(item->vert[i]->row->desc1->crd, ROTA, trans); 

        rmsd += dist_fast(c, item->vert[i]->row->desc2->crd);
      }
      
    }
    
    gsl_vector_free(vec);
    
    /* vrnemo rmsd */
    return sqrt( rmsd / unique_resi.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}

double Kabsch::superimpose_only_calpha(Item *item) {
  /*
    Funkcija prilega samo alpha atome in vrne rmsd
    med C-alpha atomi. Na ta nacin preprecis
    primere, ko sta prilegani funkcionalni skupini obrnjeni ena proti drugi.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);
  
  int st = 0;


  /* 
     dodamo unique calpha v X in Y matriko, ce je na primer 2x ARG201--ARG345, dodamo samo enkrat calpha, 
     ce je ARG201--ARG345 in ARG201--ARG346, potem dodamo dvakrat calpha.
  */
  set<pair<int, int> > unique_resi;
  unique_resi.clear();

  for (unsigned int i = 0; i < item->vert.size(); i++) {

    int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
    int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

    /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
    if (unique_resi.find( make_pair(resi1,resi2) ) == unique_resi.end() ) {

      unique_resi.insert(make_pair(resi1, resi2) );

      gsl_vector_set(vec, 0, item->vert[i]->row->desc1->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc1->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc1->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(X, st, vec);
      
      gsl_vector_set(vec, 0, item->vert[i]->row->desc2->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc2->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc2->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(Y, st, vec);

      st++;

    }
  }

  /* ce prileganje uspe, alignamo calpha atome z izracunano matriko ROTA in trans in izracunamo RMSD med calpha */
  if (kabsch(st, X, Y, ROTA, trans, NULL) == 1) {

  
    double rmsd = 0.0;

    unique_resi.clear();

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
      int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

      /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
      if (unique_resi.find( make_pair(resi1, resi2) ) == unique_resi.end() ) {

        unique_resi.insert( make_pair(resi1, resi2) );

        Coor c = rotate_vector(item->vert[i]->row->desc1->atom->acid_start->next->crd, ROTA, trans); 

        rmsd += dist_fast(c, item->vert[i]->row->desc2->atom->acid_start->next->crd);
      }
      
    }
    
    gsl_vector_free(vec);
    
    /* vrnemo rmsd */
    return sqrt( rmsd / unique_resi.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}

double Kabsch::superimpose_calpha_cbeta(Item *item) {
  /*
    Funkcija prilega alpha in beta atome in vrne rmsd
    med C-alpha atomi. Na ta nacin preprecis
    primere, ko sta prilegani funkcionalni skupini obrnjeni ena proti drugi.

    Tukaj dolocimo tudi item->U in item->t, saj ROTA kaze na item->U in trans na item->t.
  */
  gsl_matrix *ROTA =  item->U;
  gsl_vector *trans = item->t;
  gsl_vector *vec = gsl_vector_alloc(3);
  
  int st = 0;


  /* 
     dodamo unique calpha v X in Y matriko, ce je na primer 2x ARG201--ARG345, dodamo samo enkrat calpha, 
     ce je ARG201--ARG345 in ARG201--ARG346, potem dodamo dvakrat calpha.
  */
  set<pair<int, int> > unique_resi;
  unique_resi.clear();

  for (unsigned int i = 0; i < item->vert.size(); i++) {

    int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
    int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

    /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
    if (unique_resi.find( make_pair(resi1,resi2) ) == unique_resi.end() ) {

      unique_resi.insert(make_pair(resi1, resi2) );

      gsl_vector_set(vec, 0, item->vert[i]->row->desc1->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc1->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc1->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(X, st, vec);
      
      gsl_vector_set(vec, 0, item->vert[i]->row->desc2->atom->acid_start->next->crd.x);
      gsl_vector_set(vec, 1, item->vert[i]->row->desc2->atom->acid_start->next->crd.y);
      gsl_vector_set(vec, 2, item->vert[i]->row->desc2->atom->acid_start->next->crd.z);
      
      gsl_matrix_set_row(Y, st, vec);

      st++;

      /* dodajmo cbeta ce ni GLY */
      if (one_letter_code(item->vert[i]->row->desc1->atom->resn) != 'G' &&
          one_letter_code(item->vert[i]->row->desc2->atom->resn) != 'G') 
        {
          gsl_vector_set(vec, 0, item->vert[i]->row->desc1->atom->acid_start->next->next->next->next->crd.x);
          gsl_vector_set(vec, 1, item->vert[i]->row->desc1->atom->acid_start->next->next->next->next->crd.y);
          gsl_vector_set(vec, 2, item->vert[i]->row->desc1->atom->acid_start->next->next->next->next->crd.z);
      
          gsl_matrix_set_row(X, st, vec);
      
          gsl_vector_set(vec, 0, item->vert[i]->row->desc2->atom->acid_start->next->next->next->next->crd.x);
          gsl_vector_set(vec, 1, item->vert[i]->row->desc2->atom->acid_start->next->next->next->next->crd.y);
          gsl_vector_set(vec, 2, item->vert[i]->row->desc2->atom->acid_start->next->next->next->next->crd.z);
      
          gsl_matrix_set_row(Y, st, vec);

          st++;
        }

    }
  }

  /* ce prileganje uspe, alignamo calpha atome z izracunano matriko ROTA in trans in izracunamo RMSD med calpha */
  if (kabsch(st, X, Y, ROTA, trans, NULL) == 1) {

  
    double rmsd = 0.0;

    unique_resi.clear();

    for (unsigned int i = 0; i < item->vert.size(); i++) {

      int resi1 = item->vert[i]->row->desc1->atom->acid_start->next->resi;
      int resi2 = item->vert[i]->row->desc2->atom->acid_start->next->resi;

      /* upostevamo le en calpha, tudi ce vec deskriptorjev kaze nanj */
      if (unique_resi.find( make_pair(resi1, resi2) ) == unique_resi.end() ) {

        unique_resi.insert( make_pair(resi1, resi2) );

        Coor c = rotate_vector(item->vert[i]->row->desc1->atom->acid_start->next->crd, ROTA, trans); 

        rmsd += dist_fast(c, item->vert[i]->row->desc2->atom->acid_start->next->crd);
      }
      
    }
    
    gsl_vector_free(vec);
    
    /* vrnemo rmsd */
    return sqrt( rmsd / unique_resi.size() );
  }  
  else {
    gsl_vector_free(vec);
    return -1.0;
  }

}


bool Kabsch::compare_matrices(gsl_matrix *U1, gsl_matrix *U2, gsl_vector *t1, gsl_vector *t2) {

  if (fabs(gsl_matrix_get(U1, 0, 0)  - gsl_matrix_get(U2, 0, 0)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 0, 1)  - gsl_matrix_get(U2, 0, 1)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 0, 2)  - gsl_matrix_get(U2, 0, 2)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 1, 0)  - gsl_matrix_get(U2, 1, 0)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 1, 1)  - gsl_matrix_get(U2, 1, 1)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 1, 2)  - gsl_matrix_get(U2, 1, 2)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 2, 0)  - gsl_matrix_get(U2, 2, 0)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 2, 1)  - gsl_matrix_get(U2, 2, 1)  ) < MATRIX_CUTOFF &&
      fabs(gsl_matrix_get(U1, 2, 2)  - gsl_matrix_get(U2, 2, 2)  ) < MATRIX_CUTOFF &&

      fabs(gsl_vector_get(t1, 0)     - gsl_vector_get(t2, 0)     ) < VECTOR_CUTOFF &&
      fabs(gsl_vector_get(t1, 1)     - gsl_vector_get(t2, 1)     ) < VECTOR_CUTOFF &&
      fabs(gsl_vector_get(t1, 2)     - gsl_vector_get(t2, 2)     ) < VECTOR_CUTOFF)
    return true;


  return false;

}



//    double Kabsch::rmsd(Item *item, gsl_matrix *ROTA, gsl_vector *trans) {
//    
//      double rez=0;
//      gsl_vector *y = gsl_vector_alloc(3);
//      gsl_vector *y1 = gsl_vector_alloc(3);
//      gsl_vector *x = gsl_vector_alloc(3);
//    
//      for (unsigned int j = 0; j < item->vert.size(); j++) {
//        gsl_matrix_get_row(x, X, j);
//    
//        gsl_blas_dgemv(CblasNoTrans, 1, ROTA, x, 0, y);
//    
//        gsl_vector_add(y, trans);   // ali je to manjkalo ?
//    
//        gsl_matrix_get_row(y1, Y, j);
//        gsl_vector_sub(y, y1);
//        gsl_vector_mul(y, y);
//        rez += sqrt(gsl_vector_get(y, 0) + gsl_vector_get(y, 1) + gsl_vector_get(y, 2));
//      }
//      gsl_vector_free(y);
//      gsl_vector_free(y1);
//      gsl_vector_free(x);
//      return rez / item->vert.size();
//    }

Coor Kabsch::rotate_vector(Coor crd, gsl_matrix *ROTA, gsl_vector *trans) {
  /*
    Rotate a vector in crd according to ROTA and trans and return the rotated vector.
  */
  Coor c;
  gsl_vector *vec = gsl_vector_alloc(3);
  gsl_vector *y = gsl_vector_alloc(3);
  
  gsl_vector_set(vec, 0, crd.x);
  gsl_vector_set(vec, 1, crd.y);
  gsl_vector_set(vec, 2, crd.z);
  
//  gsl_vector_add(vec, trans);
//
//  gsl_blas_dgemv(CblasTrans, 1, ROTA, vec, 0, y);
  gsl_blas_dgemv(CblasNoTrans, 1, ROTA, vec, 0, y);
  
  gsl_vector_add(y, trans);
  
  c.x = gsl_vector_get(y, 0);
  c.y = gsl_vector_get(y, 1);
  c.z = gsl_vector_get(y, 2);
  
  gsl_vector_free(vec);
  gsl_vector_free(y);
  return c;
}


Coor Kabsch::rotate_vector_inv(Coor crd, gsl_matrix *ROTA, gsl_vector *trans) {
  /*
    Inversely rotate a vector in crd according to ROTA and trans and return the rotated vector.
    The problem is that in mysql are stored rotation and translation of PROTEIN1 to PROTEIN2. To rotate
    PROTEIN2 to PROTEIN1, the equation has to be inverted to get vec: ROTA*vec+trans = y | * A^(-1)
    A^(-1)=A(transposed) , therefore to get back the initial vector, vec = ROTA(transposed) * (y - trans)
  */
  Coor c;
  gsl_vector *vec = gsl_vector_alloc(3);
  gsl_vector *y = gsl_vector_alloc(3);

  gsl_vector_set(vec, 0, crd.x);
  gsl_vector_set(vec, 1, crd.y);
  gsl_vector_set(vec, 2, crd.z);
  
  gsl_vector_sub(vec, trans);

  gsl_blas_dgemv(CblasTrans, 1, ROTA, vec, 0, y);

  c.x = gsl_vector_get(y, 0);
  c.y = gsl_vector_get(y, 1);
  c.z = gsl_vector_get(y, 2);

  gsl_vector_free(vec);
  gsl_vector_free(y);
  return c;
}

void Kabsch::invert(gsl_matrix *ROTA, gsl_vector *trans) {
  /*
    Transponiramo matriko in negiramo vektor (in-place).
    y = ROTA*x + trans --> x = ROTAt*y - ROTAt*trans
    ROTA --> ROTAt
    trans --> ROTAt * (-trans)
    
  */

  gsl_matrix_transpose(ROTA);

  gsl_vector *vec = gsl_vector_alloc(3);
  gsl_vector_set_all(vec, 0.0);

  gsl_vector_sub(vec, trans);
  gsl_blas_dgemv( CblasNoTrans, 1, ROTA, vec, 0.0, trans );  

  gsl_vector_free(vec);
}

void Kabsch::output(gsl_matrix *ROTA, gsl_vector *trans) {
  /*
    Izpisemo vektor in matriko.
  */
  cout << "Vektor = (" << gsl_vector_get(trans, 0) << "," << gsl_vector_get(trans, 1) << "," << gsl_vector_get(trans, 2) << ")" << endl;
  cout << "Matrika = (" << gsl_matrix_get(ROTA, 0, 0) << "," << gsl_matrix_get(ROTA, 0, 1) << "," << gsl_matrix_get(ROTA, 0, 2) << ")" << endl;
  cout << "          (" << gsl_matrix_get(ROTA, 1, 0) << "," << gsl_matrix_get(ROTA, 1, 1) << "," << gsl_matrix_get(ROTA, 1, 2) << ")" << endl;
  cout << "          (" << gsl_matrix_get(ROTA, 2, 0) << "," << gsl_matrix_get(ROTA, 2, 1) << "," << gsl_matrix_get(ROTA, 2, 2) << ")" << endl;
}

void Kabsch::set_rota_trans(gsl_matrix *ROTA, gsl_vector *trans, char **row) {
  /*
    Nastavimo rotacijsko matriko in vektor iz char arraya.
  */
    gsl_matrix_set(ROTA, 0, 0, atof(row[0]) );
    gsl_matrix_set(ROTA, 0, 1, atof(row[1]) );
    gsl_matrix_set(ROTA, 0, 2, atof(row[2]) );
    gsl_matrix_set(ROTA, 1, 0, atof(row[3]) );
    gsl_matrix_set(ROTA, 1, 1, atof(row[4]) );
    gsl_matrix_set(ROTA, 1, 2, atof(row[5]) );
    gsl_matrix_set(ROTA, 2, 0, atof(row[6]) );
    gsl_matrix_set(ROTA, 2, 1, atof(row[7]) );
    gsl_matrix_set(ROTA, 2, 2, atof(row[8]) );

    gsl_vector_set(trans, 0, atof(row[9]) );
    gsl_vector_set(trans, 1, atof(row[10]) );
    gsl_vector_set(trans, 2, atof(row[11]) );
}

//void Kabsch::set_rota_trans(gsl_matrix *ROTA, gsl_vector *trans, vector<string> row) {
//  /*
//    Nastavimo rotacijsko matriko in vektor iz char arraya.
//  */
//    gsl_matrix_set(ROTA, 0, 0, atof(row[0].c_str()) );
//    gsl_matrix_set(ROTA, 0, 1, atof(row[1].c_str()) );
//    gsl_matrix_set(ROTA, 0, 2, atof(row[2].c_str()) );
//    gsl_matrix_set(ROTA, 1, 0, atof(row[3].c_str()) );
//    gsl_matrix_set(ROTA, 1, 1, atof(row[4].c_str()) );
//    gsl_matrix_set(ROTA, 1, 2, atof(row[5].c_str()) );
//    gsl_matrix_set(ROTA, 2, 0, atof(row[6].c_str()) );
//    gsl_matrix_set(ROTA, 2, 1, atof(row[7].c_str()) );
//    gsl_matrix_set(ROTA, 2, 2, atof(row[8].c_str()) );
//
//    gsl_vector_set(trans, 0, atof(row[9].c_str()) );
//    gsl_vector_set(trans, 1, atof(row[10].c_str()) );
//    gsl_vector_set(trans, 2, atof(row[11].c_str()) );
//}

