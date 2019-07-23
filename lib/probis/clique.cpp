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

#include "clique.h"
#include "debug.hpp"
#include "item.h"
#include "subgraph.h"

void Clique::clear_connect(Item* item) {
    /*
      Clears the connect table e[][] of edges of m-th subgraph .
    */

    for (unsigned int i = 0; i < item->edge.size(); i++)
        set_zero(item->edge[i].first, item->edge[i].second);
}

void Clique::wipe_connect() {
    for (int i = 0; i < MAX_VERTICES / WORD; i++)
        for (int j = 0; j < MAX_VERTICES; j++) set_zero_fast(i, j);
}

void Clique::assign_connect(Item* item) {
    /*
      Assign edges of i-th subgraph to the connect table e[][]; when an edge is
      present between
      vertices i and j, e[i][j] is set to 1, otherwise 0.
    */

    for (unsigned int i = 0; i < item->edge.size(); i++)
        set_one(item->edge[i].first, item->edge[i].second);

    //  for (int j = 0; j < item->size; j++)
    //    for (int k = 0; k < item->vert[j].item.size; k++)
    //      set_one(item->vert[j].vertex, item->vert[j].item.vert[k].vertex);
}

bool by_degree(Vert* vi, Vert* vj) {
    return (vi->size_edge > vj->size_edge);
}  // callback for sorting vertices by degree

// void Clique::sort_vertices(Item item) {
//  /*
//    Sort vertices V in a descending order with respect to their degree
//  */
//  int k;
//  Vert tmp;
//  for (int i = 0; i < item.size; i++) {
//    k = i;
//    for (int j = i + 1; j < item.size; j++)
//      if (item.vert[k].item.size < item.vert[j].item.size) {
//        k = j;
//      }
//    tmp = item.vert[i];
//    item.vert[i] = item.vert[k];
//    item.vert[k] = tmp;
//  }
//}

// int Clique::CUT1(int p, Vert2 *B) {
//  /*
//    return 1 if intersection of A and B is not empty
//    return 0 if there are no elements in common to A and B
//    in C we return the intersecting elements
//  */
//  Vert2 *tmp = B;
//  for (; tmp != NULL; tmp = tmp->next)
//    if (get(p, tmp->vert.vertex))
//      break;
//  if (tmp == NULL)
//    return 0;
//  else
//    return 1;
//}

int Clique::CUT1(int p, vector<Vert*>& B) {
    /*
      return 1 if intersection of A and B is not empty
      return 0 if there are no elements in common to A and B
      in C we return the intersecting elements
    */
    unsigned int j;

    for (j = 0; j < B.size(); j++) {
#ifndef NDEBUG
        dbgmsg("VERB> CUT1 B[" << j << "] = " << B[j]->vertex);
#endif
        if (get(p, B[j]->vertex)) break;
    }
    if (j == B.size())
        return 0;
    else
        return 1;
}

int Clique::CUT2(int p, vector<Vert*>& B, vector<Vert*>& A) {
    /*
      return 1 if intersection of A and B is not empty
      return 0 if there are no elements in common to A and B
      in C we return the intersecting elements
    */
    //  unsigned int j;
    //  A.size = 0;
    for (int j = 0; j < (signed)B.size() - 1; j++) {
        if (get(p, B[j]->vertex)) {
            //      A.push_back( B[j] ); // (vi, color, edges, row)
            A.push_back(new Vert(B[j]->lista_id, B[j]->cluster_id, B[j]->vertex,
                                 0, B[j]->size_edge,
                                 B[j]->row));  // (vi, color, edges, row)
            //      A.vert[A.size++] = B.vert[j];
            //      RESULT.vert[RESULT.size].row = B.vert[j].row;
            //      RESULT.vert[RESULT.size++].vertex = B.vert[j].vertex;
        }
    }

    if (A.size() == 0)
        //  if (RESULT.size == 0)
        return 0;
    else
        return 1;
}

void Clique::COLOR_SORT(vector<Vert*>& R) {
    unsigned int i;
    int j, k, p, min_k, maxno;
    vector<vector<Vert*> > C;
    C.reserve(R.size() + 2);
#ifndef NDEBUG
    dbgmsg("VERB> COLOR_SORT C.size() = " << C.size());
#endif

    maxno = 1;
    min_k = QMAX->vert.size() - Q->vert.size();  // change to + 1 when
                                                 // calculating only one maximum
                                                 // clique !!

#ifndef NDEBUG
    dbgmsg("VERB> COLOR_SORT min_k = " << min_k);
    dbgmsg("VERB> COLOR_SORT R.size() = " << R.size());
#endif
    vector<Vert*> vec;
    C.push_back(vec);
    C.push_back(vec);
    C.push_back(vec);
    C[1].reserve(R.size());
    C[2].reserve(R.size());  // for sure reserve a positive number (not
                             // R.vert.size() - 1 )

    i = 0;
    j = 0;
    while (i < R.size()) {
        p = R[i]->vertex;
        k = 1;
#ifndef NDEBUG
        dbgmsg("VERB> COLOR_SORT C[" << k << "].size() = " << C[k].size());
#endif
        while (CUT1(p, C[k]) != 0) k++;

        if (k > maxno) {
            maxno = k;
            C.push_back(vec);
            C[maxno + 1].reserve(R.size());
        }

        C[k].push_back(R[i]);

        if (k < min_k) {
            R[j++] = R[i];
        }
        i++;
    }
    if (j > 0) R[j - 1]->color = 0;
    if (min_k <= 0) min_k = 1;
    for (k = min_k; k <= maxno; k++)
        for (i = 0; i < C[k].size(); i++) {
            R[j] = C[k][i];
            R[j++]->color = k;
        }
    C.clear();
}

// void Clique::NUMBER_SORT(Item R) {
//  int i, j, k;
//  int min_k, p;
//  vector<vector< Vert > > C ;
//
////  Vert2 *L = NULL, *c, *tmp = NULL, *min_c = NULL;
//
//  min_k = QMAX->vert.size() - Q.vert.size() + 1;
//  i = 0;
//  j = 0;
//  while (i < R.vert.size()) {
//    p = R.vert[i]->vertex;
//    k = 1;
//    //    c = L;
//    while (CUT1(p, C[k] ) != 0) {
//      k++;
//      tmp = c;
//      c = c->next_level;
//    }
//    if (c == NULL) {
//      if (i == 0) {
//        L = new Vert2();
//        L->vert = R.vert[i];
//        L->last = L;
//        c = L;
//      }
//      else {
//        tmp->next_level = new Vert2();
//        tmp->next_level->vert = R.vert[i];
//        tmp->next_level->last = tmp->next_level;
//        c = tmp->next_level;
//      }
//    }
//    else {
//      c->last->next = new Vert2();
//      c->last->next->vert = R.vert[i];
//      c->last = c->last->next;
//    }
//    if (k < min_k) {
//      R.vert[j++] = R.vert[i];
//    }
//    if (k == min_k)
//      min_c = c;
//
//    i++;
//  }
//  if (j > 0) R.vert[j-1].color = 0;
//  if (min_k <= 0) { min_c = L; min_k = 1; }
//  k = min_k;
//  for (c = min_c; c != NULL; c = c->next_level) {
//    for (tmp = c; tmp != NULL; tmp = tmp->next) {
//      R.vert[j] = tmp->vert;
//      R.vert[j++].color = k;
//    }
//    k++;
//  }
//  L->free_vert2();
//}

void Clique::COPY(vector<Vert*>& A, vector<Vert*>& B) {
    /*
      Deallocate each vertex of A and clear the vector array, then refill it
      with copies of vertices of B
    */
    for (unsigned int i = 0; i < A.size(); i++) delete A[i];
    A.clear();
    A.reserve(B.size());
    for (unsigned int i = 0; i < B.size(); i++)
        A.push_back(new Vert(B[i]->lista_id, B[i]->cluster_id, B[i]->vertex,
                             B[i]->color, B[i]->size_edge, B[i]->row));
}

void Clique::DELETE(vector<Vert*>& A) {
    /*
      Deallocate each vertex of A and clear the vector array
    */
    for (unsigned int i = 0; i < A.size(); i++) delete A[i];
    A.clear();
}

// void Clique::DEGREES_SORT(vector<Vert*> &R) {
//
//  for (unsigned int i = 0; i < R.size(); i++) {
//    R[i]->size_edge = 0;
//    for (unsigned int j = 0; j < i; j++)
//      if (get(R[i]->vertex, R[j]->vertex)) {
//        R[i]->size_edge++;
//        R[j]->size_edge++;
//      }
//  }
//
//  sort(R.begin(), R.end(), by_degree );
//}

//  void Clique::set_size_edge(vector<Vert* > &R) {
//
//    for (unsigned int i = 0; i < R.size(); i++) {
//      R[i]->size_edge = 0;
//      for (unsigned int j = 0; j < i; j++)
//        if (get(R[i]->vertex, R[j]->vertex)) {
//          R[i]->size_edge++;
//  //        R.vert[j]->size_edge++;
//        }
//    }
//
//  //  sort(R.vert.begin(), R.vert.end(), by_degree );
//  //  sort_vertices(R);
//  }

// int Clique::weight(Item *A) {
//#ifndef NDEBUG
//  dbgmsg("Clique::weight A->vert.size() = " << A->vert.size());
//#endif
//  int w = 0;
//  for (int i = 0; i < (signed) A->vert.size() - 1; i++) {   // signed here is
//  crucial otherwise expression is INF
//    for (unsigned int j = i; j < A->vert.size(); j++)
////      w += lookup_weight[pair<int, int > (make_pair(A->vert[i]->vertex,
///A->vert[j]->vertex ) ) ];
//      w += A->weight_edge[pair<int, int > (make_pair(A->vert[i]->vertex,
//      A->vert[j]->vertex ) ) ];
//  }
//#ifndef NDEBUG
//  dbgmsg("Clique::weight w = " << w);
//#endif
//  return w;
//}

void Clique::EXPAND(vector<Vert*>& R, int level) {
    /*
      First level is level = 1;
    */
    int p;
    vector<Vert*> Rp;

    //  double tmp = 0.0;

    //  CNT.vert[level]->vertex = CNT.vert[level]->vertex + CNT.vert[level -
    //  1]->vertex - CNT.vert[level]->color;
    //  CNT.vert[level]->color = CNT.vert[level - 1]->vertex;

    while (R.size() != 0) {
        p = R.back()->vertex;  // p is assigned the last vertex in R
        //    p = R.vert[R.vert.size() - 1]->vertex;  // p is assigned the last
        //    vertex in R

        if (Q->vert.size() + R.back()->color >= QMAX->vert.size()) {
            //    if (Q->vert.size() + R.back()->color > QMAX->vert.size() &&
            //    Q->vert.size() < 5) {

            //    if (Q->vert.size() + R.vert[R.vert.size() - 1]->color >
            //    QMAX->vert.size() ) {

            //      Q->vert.push_back( R.back() );
            Q->vert.push_back(new Vert(R.back()->lista_id, R.back()->cluster_id,
                                       p, 0, 0, R.back()->row));

// Q->vert.push_back(new Vert(p, 0, 0, R.vert[R.vert.size() - 1]->row ) );

#ifndef NDEBUG
            dbgmsg("VERB> EXPAND p = " << p);
#endif
            //      exit(1);
            //      Q.vert[Q.size].row = R.vert[R.size - 1].row;
            //      Q.vert[Q.size++].vertex = p;

            //      Q.rmsd_old = Q.rmsd;

            //      gsl_matrix_memcpy(Q.U_old, Q.U);
            //      gsl_vector_memcpy(Q.t_old, Q.t);

            //      if ((tmp = superimpose(Q, Q.U, Q.t)) < RMSD_INCR + 1000.0) {
            //      // ocitno je preverjanje rmsd-jev izklopljeno!
            //
            //        Q.rmsd = tmp;

            Rp.reserve(R.size());
            //        Rp.vert = (Vert*) calloc(R.size, sizeof(Vert));

            if (CUT2(p, R, Rp) != 0) {
//          if ((double)CNT.vert[level]->vertex/pk < num_level) {
//  //          nk++;
//  //          nk_size+=Rp.size();
//            DEGREES_SORT(Rp);
//  //#ifndef NDEBUG
//            dbgmsg("VERB> EXPAND za degrees_sort");
//  //#endif
//          }
#ifndef NDEBUG
                dbgmsg("VERB> level = " << level);
                dbgmsg("VERB> EXPAND R.size() = " << R.size());
                for (unsigned int j = 0; j < R.size(); j++) {
                    cout << R[j]->vertex << " (" << R[j]->size_edge << ")";
                    //          for (unsigned int k = 0; k < R.edge.size(); k++)
                    //          {
                    //            if (R.edge[k].first == R.vert[j]->vertex )
                    //              cout << R.edge[k].second << " ";
                    //            //      if (item->edge[k].second ==
                    //            item->vert[j]->vertex )
                    //            //        cout << item->edge[k].first << " ";
                    //          }
                    cout << endl;
                }
                dbgmsg("VERB> EXPAND Rp.size() = " << Rp.size());
                for (unsigned int j = 0; j < Rp.size(); j++) {
                    cout << Rp[j]->vertex << " (" << Rp[j]->size_edge << ")";
                    //          for (unsigned int k = 0; k < Rp.edge.size();
                    //          k++) {
                    //            if (Rp.edge[k].first == Rp.vert[j]->vertex )
                    //              cout << Rp.edge[k].second << " ";
                    //            //      if (item->edge[k].second ==
                    //            item->vert[j]->vertex )
                    //            //        cout << item->edge[k].first << " ";
                    //          }
                    cout << endl;
                }
#endif

                COLOR_SORT(Rp);
                //        exit(1);

                //        CNT.vert[level]->vertex++;
                pk++;
                /* na vsakih 100 korakov preveri, ce ni potekel cas */
                if (!(pk % 100))
                    if ((double)(clock() - start) / CLOCKS_PER_SEC >
                        CLIQUE_TLIMIT) {
                        dbgmsg("Clique::EXPAND Time limit breached after "
                               << (double)(clock() - start) / CLOCKS_PER_SEC
                               << "s");

                        end_condition = true;
                    }

                if (!end_condition) EXPAND(Rp, level + 1);
            }
            //      else if (Q->vert.size() == 2 ) {
            else if (Q->vert.size() >= QMAX->vert.size()) {
                //        /* if a larger clique was found than one in QMAX
                //        delete all QMAX */
                //        if (Q->vert.size() > QMAX.back()->vert.size() ) {
                //          for (unsigned int i = 0; i < QMAX.size(); i++) //
                //          delete all QMAX vector
                //            delete QMAX[i];
                //          QMAX.clear();
                //        }
                //        //        QMAX->vert = Q->vert;
                //        QMAX.push_back(new Item);
                //        if (Q->vert.size() == 2 ) {

                if (Q->vert.size() == QMAX->vert.size()) {
                    //          for (unsigned int i = 0; i < Q->vert.size();
                    //          i++)
                    //            cout << Q->vert[i]->vertex << " " ;
                    //          cout << endl;

                    //          if ( _state == RESULTS ) {
                    //            if ( weight(Q) > weight(QMAX) )
                    //              COPY(QMAX->vert, Q->vert);
                    //          }
                    //          else
                    //            if ( superimpose(*Q, Q->U, Q->t) <
                    //            superimpose(*QMAX, QMAX->U, QMAX->t) ) {
                    if (superimpose_calpha(Q) < superimpose_calpha(QMAX)) {
#ifndef NDEBUG
                        dbgmsg("Q->rmsd = " << superimpose_calpha(Q));
                        dbgmsg("QMAX->rmsd = " << superimpose_calpha(QMAX));
#endif
                        //              dbgmsg("Q->rmsd = " << superimpose(*Q,
                        //              Q->U, Q->t));
                        //              dbgmsg("QMAX->rmsd = " <<
                        //              superimpose(*QMAX, QMAX->U, QMAX->t));
                        COPY(QMAX->vert, Q->vert);
                    }
                } else {
                    COPY(QMAX->vert, Q->vert);
                }
#ifndef NDEBUG
                dbgmsg("EXPAND: Q->vert.size() = " << Q->vert.size()
                                                   << " QMAX->vert.size() = "
                                                   << QMAX->vert.size());
#endif

#ifndef NDEBUG
                dbgmsg("found maximal clique QMAX->vert.size() = "
                       << QMAX->vert.size());
                cout << "QMAX elements are: ";
                for (unsigned int i = 0; i < QMAX->vert.size(); i++)
                    cout << QMAX->vert[i]->vertex << " ";
                cout << endl;
#endif
                //        QMAX.push_back(new Item() );
            }
            DELETE(Rp);
            //      /* Rp will be deallocated automatically when leaving EXPAND
            //      */
            //      for (int i = 0; i < Rp.size(); i++ )
            //        delete Rp[i];
            //      Rp.clear();
            //      delete [] Rp.vert;

            //      }
            //      else {
            //        // in Q.rmsd is old rmsd without vertex p
            //        if (Q.size - 1 >= QMAX->size) {
            //          if (Q.size - 1 > QMAX->size) {
            //            Q.size--;
            //#ifndef NDEBUG
            //            dbgmsg("morebitna napaka Q.rmsd = " << Q.rmsd << " tmp
            //            = " << tmp << " super = " << superimpose(Q, Q.U,
            //            Q.t));
            //#endif
            //            Q.size++;
            //            gsl_matrix_memcpy(Q.U, Q.U_old);
            //            gsl_vector_memcpy(Q.t, Q.t_old);
            //            COPY(QMAX, Q);
            //            QMAX->size = Q.size - 1;
            //          }
            //          else if (Q.rmsd < QMAX->rmsd) {
            //            Q.size--;
            //#ifndef NDEBUG
            //            dbgmsg("morebitna napaka 2 Q.rmsd = " << Q.rmsd << "
            //            tmp = " << tmp << " super = " << superimpose(Q, Q.U,
            //            Q.t));
            //#endif
            //            Q.size++;
            //            gsl_matrix_memcpy(Q.U, Q.U_old);
            //            gsl_vector_memcpy(Q.t, Q.t_old);
            //            COPY(QMAX, Q);
            //            QMAX->size = Q.size - 1;
            //          }
            //        }
            //      }
            delete Q->vert.back();
            Q->vert.pop_back();
            //      Q.size--;
        } else {
            return;
        }
        delete R.back();
        R.pop_back();
        //    R.size--;
    }
}

void Clique::MCQ(Item* item, Item* qmax) {
    int max_degree;

    QMAX = qmax;

    Q = new Item();

    //  Q->weight_edge = QMAX->weight_edge;

    pk = 0;

    end_condition = false;

    Q->vert.reserve(item->vert.size());

    sort(item->vert.begin(), item->vert.end(), by_degree);

    assign_connect(item);

#ifndef NDEBUG
    for (unsigned int j = 0; j < item->vert.size(); j++) {
        cout << item->vert[j]->vertex << " (" << item->vert[j]->size_edge
             << ") : ";
        for (unsigned int k = 0; k < item->edge.size(); k++) {
            if (item->edge[k].first == item->vert[j]->vertex)
                cout << item->edge[k].second << " ";
        }
        cout << endl;
    }
#endif

    max_degree = item->vert[0]->size_edge;

#ifndef NDEBUG
    dbgmsg("max degree = " << max_degree);
#endif

    for (int j = 0; j < max_degree; j++) item->vert[j]->color = j + 1;
    for (unsigned int j = max_degree; j < item->vert.size(); j++)
        item->vert[j]->color = max_degree + 1;

    vector<Vert*> R;

    COPY(R, item->vert);

    EXPAND(R, 1);

#ifndef NDEBUG
    dbgmsg("QMAXFINAL->rmsd = " << superimpose_calpha(QMAX));
#endif

    DELETE(R);  // delete vertices and deallocate each Vert* !! )

    QMAX->calpha_rmsd = superimpose_calpha(QMAX);

#ifndef NDEBUG
    dbgmsg("QMAX->vert.size() = " << QMAX->vert.size());
    //~ dbgmsg("QMAX->rmsd = " << QMAX->rmsd);
    //  dbgmsg("QMAX->weight = " << QMAX->weight);
    dbgmsg("Number of steps = " << pk);
    cout << "QMAX elements are: ";
    for (unsigned int i = 0; i < QMAX->vert.size(); i++)
        cout << QMAX->vert[i]->vertex << " ";
    cout << endl;
#endif
    //  out->output_align4(QMAX);
    //  out->output_align_protein(mol1, mol2, QMAX);
    //  qmax = QMAX.back();

    clear_connect(item);

    //  CNT.vert.clear();

    //  Q.vert.clear();
    delete Q;  // item destructor destroys vert

    //  for (unsigned int i = 0; i < QMAX.size(); i++) // delete all QMAX vector
    //    delete QMAX[i];
    //  QMAX.clear();

    //

    //  delete [] CNT.vert;
    //  delete [] Q.vert;

    //  gsl_matrix_free(Q.U);
    //  gsl_matrix_free(Q.U_old);
    //  gsl_vector_free(Q.t);
    //  gsl_vector_free(Q.t_old);
}

void Clique::max_clique(Subgraph* subgraph) {
    dbgmsg("CLIQUE> Cliques = " << subgraph->item.size());
#ifndef NDEBUG
    clock_t lokalni_start = clock();
#endif

    for (unsigned int i = 0; i < subgraph->item.size(); i++) {
#ifndef NDEBUG
        dbgmsg(i << "-th clique");
        dbgmsg(i << "-th subgraph size = " << subgraph->item[i]->vert.size());
#endif
        start = clock();  // ne spreminjaj, mora bit tukaj - glej EXPAND !

        MCQ(subgraph->item[i], subgraph->qmax[i]);
#ifndef NDEBUG
        dbgmsg("CLIQUE> Finished "
               << i + 1 << "/" << subgraph->item.size() << " cliques ["
               << (double)(clock() - start) / CLOCKS_PER_SEC << "s]");
#endif
    }

    dbgmsg("Clique::max_clique() Time = " << (double)(clock() - lokalni_start) /
                                                 CLOCKS_PER_SEC);
}

void Clique::read_dimacs(Subgraph* subgraph, char name[], int num_i) {
    ifstream f(name);
    char buffer[256], token[20];
    unsigned int i, j;
    int vi, vj;

    if (!f.is_open()) {
        throw Err("Error (CLIQUE) : Cannot open dimacs format file!", 22);
    }

    wipe_connect();

    //  subgraph->size = num_i + 1;
    subgraph->item.reserve(num_i + 1);
    //  subgraph->item = (Item*) calloc(num_i + 1, sizeof(Item));

    int num_edges = 0;

    while (!f.eof()) {
        f.getline(buffer, 250);
        if (buffer[0] == 'p') {
            dbgmsg(buffer);
            unsigned int size_vertex;
            sscanf(&buffer[7], "%d", &size_vertex);
            dbgmsg("Num of vertices in graph = " << size_vertex);
            subgraph->item[num_i]->vert.reserve(size_vertex);
            subgraph->item[num_i]->edge.reserve(
                3 * size_vertex);  // we reserve a bit more space for edges (but
                                   // not exactly)

            //      subgraph->item[num_i].vert = (Vert*)
            //      calloc(subgraph->item[num_i].size, sizeof(Vert));
            for (i = 0; i < size_vertex; i++) {
                subgraph->item[num_i]->vert.push_back(new Vert(i, 0, 0, NULL));
                //        subgraph->item[num_i].vert[i].vertex = i;
                //        subgraph->item[num_i].vert[i].item.size = 0;
                //        subgraph->item[num_i].vert[i].item.vert = (Vert*)
                //        calloc(subgraph->item[num_i].size, sizeof(Vert));
            }
        }
        if (buffer[0] == 'e') {
            i = 2;
            j = 0;
            while (buffer[i] != ' ') {
                token[j++] = buffer[i];
                i++;
            }
            token[j] = '\0';
            vi = atoi(token);
            i++;
            j = 0;
            while (buffer[i] != ' ') {
                token[j++] = buffer[i];
                i++;
            }
            token[j] = '\0';
            vj = atoi(token);
            vi--;
            vj--;
            subgraph->item[num_i]->edge.push_back(make_pair(vi, vj));
            subgraph->item[num_i]->edge.push_back(make_pair(vj, vi));
            //      subgraph->item[num_i].vert[vi].item.vert[subgraph->item[num_i].vert[vi].item.size++].vertex
            //      = vj;
            //      subgraph->item[num_i].vert[vj].item.vert[subgraph->item[num_i].vert[vj].item.size++].vertex
            //      = vi;

            num_edges++;
        }
    }
    assign_connect(subgraph->item[num_i]);
#ifndef NDEBUG
    double dens =
        (double)num_edges / (subgraph->item[num_i]->vert.size() *
                             (subgraph->item[num_i]->vert.size() - 1) / 2);
    dbgmsg("|E| = " << num_edges
                    << "  |V| = " << subgraph->item[num_i]->vert.size()
                    << "  gostota = " << dens);
#endif
    f.close();
}
