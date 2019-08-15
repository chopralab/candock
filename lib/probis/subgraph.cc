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

#include "subgraph.h"
#include "item.h"
#include "vert.h"
#include "molecule.h"
#include "atom.h"
#include "product.h"
#include "kabsch.h"
#include "output.h"
#include "desc.h"
#include "score.h"
#include "debug.hpp"

/* callback functions cannot be class members */

bool by_hssp_length (const pair< TwoChResi, int> &i, const pair< TwoChResi, int> &j) { return (i.second > j.second); }

bool by_size (Item *i, Item *j) { return (i->vert.size() > j->vert.size()); }

bool by_cluster_score (Item *i, Item *j) { return (i->cluster_score > j->cluster_score); }
bool by_capacity (Item *i, Item *j) { return (i->vert.capacity() > j->vert.capacity()); }


bool by_size_blosum (Item *i, Item *j) { return ((i->vert.size() == j->vert.size() ) && (i->blosum_score < j->blosum_score ) ); }
bool by_blosum (Item *i, Item *j) { return (i->blosum_score < j->blosum_score ) ; }

bool by_size_surf_vect (Item *i, Item *j) { return (i->vert.size() == j->vert.size() && i->surf_vector_angle < j->surf_vector_angle); }

bool by_vertex (Vert *i, Vert *j) { return (i->vertex < j->vertex); }

/* end callback */


void Subgraph::extend_all(Score *score, Molecule *mol1, Molecule *mol2, Product *product, Descriptor*& desc1, Descriptor*& desc2) {
  /*
    Gremo po vseh klastrih in vsakega poskusimo podaljsati!
  */
#ifndef NDEBUG
  clock_t start = clock();
#endif
  for (int clNum = 0; clNum < (signed int) clus.size(); clNum++) {
    extend_alignments(score, mol1, mol2, product, desc1, desc2, clNum);
  }
  
  dbgmsg("Subgraph::extend_all() Time = " << (double) (clock() - start) / CLOCKS_PER_SEC );
  
}

void Subgraph::extend_first_some(Score *score, Molecule *mol1, Molecule *mol2, Product *product, Descriptor*& desc1, Descriptor*& desc2) {
  /*
    Gremo po prvih 10 najvecjih klastrih in vsakega poskusimo podaljsati! 
    Klastri morajo biti pred tem sortirani po velikosti!
  */
#ifndef NDEBUG
  clock_t start = clock();
#endif
  for (int clNum = 0; clNum < (signed int) clus.size() && clNum < 10; clNum++) {
    extend_alignments(score, mol1, mol2, product, desc1, desc2, clNum);
  }
  
  dbgmsg("Subgraph::extend_all() Time = " << (double) (clock() - start) / CLOCKS_PER_SEC );
  
}

void Subgraph::extend_alignments(Score *score, Molecule *m1, Molecule *m2, Product *p, Descriptor*& desc1, Descriptor*& desc2, int clNum) {
  /*
    Vsakemu klastru dodamo pare prileganih aminokislin, ki so slabse strukturno ohranjene.
  */
//  dbgmsg("test" );

  Kabsch *kabsch = new Kabsch();


  map <ChResi, Atom*> ca2;   // <chain_id, resi> --> atm

  /* najdemo c-alpha atome molekule 2 */
  for (Atom *atm = m2->atom; atm != NULL; atm = atm->next) {
//    dbgmsg("test = " << atm->tag );

      if (!strcmp(atm->tag, "CA") ) {

#ifndef NDEBUG
        dbgmsg("Subgraph::extend_alignment Koordinate m2 " << atm->crd.x << "," << atm->crd.y << "," << atm->crd.z );
#endif


        ca2[ ChResi(atm->chain_id, atm->resi) ] = atm;
      }
  }

  map<TwoChResi, Row*> row_map;

  /* 1. upostevamo cluster stevilka clNum */
  map <ChResi, Atom*> ca1;
  
  ca1.clear();
  
  /* 2. rotiramo m1 na m2, tako kot narekuje ROTA in TRANS v clus[clNum]  */
  for (Atom *atm = m1->atom; atm != NULL; atm = atm->next) {
    
    /* 3. upostevamo le backbone C-alpha (CA) atome obeh proteinov */
    if (!strcmp(atm->tag, "CA") ) {

#ifndef NDEBUG
      dbgmsg("Subgraph::extend_alignment Koordinate pred rotacijo " << atm->crd.x << "," << atm->crd.y << "," << atm->crd.z );
#endif
      Coor crot = kabsch->rotate_vector(atm->crd, clus[clNum]->U, clus[clNum]->t);
        
      atm->crd = crot;
      
#ifndef NDEBUG
      dbgmsg("Subgraph::extend_alignment Koordinate po rotaciji " << atm->crd.x << "," << atm->crd.y << "," << atm->crd.z );
#endif

      ca1[ ChResi(atm->chain_id,atm->resi) ] = atm;
    }
  }        
  /* HSSP je high scoring sequence pair */
  vector < pair < TwoChResi, int > > HSSP;  // (prvi resi v m1, prvi resi v m2 in dolzina)

  HSSP.clear();

  /* 4. gremo po tri aminokisline naenkrat na m1 premikamo se od zacetka proti koncu verige */
  for (map<ChResi, Atom*>::iterator it = ca1.begin(); it != --(--ca1.end()); it++ ) {

    Atom *at10, *at11, *at12;
    
    map<ChResi, Atom*>::iterator hit = it;

    at10 = hit->second;
    
    hit++; at11 = hit->second;
    hit++; at12 = hit->second;
    
    /* tri zaporedne aminokisline in vse pripadajo isti verigi */
    if ((at10->chain_id == at11->chain_id) && (at11->chain_id == at12->chain_id) && (at10->resi == at11->resi - 1) && (at11->resi == at12->resi - 1)) {


#ifndef NDEBUG
      dbgmsg("Subgraph::extend_alignment Tri zaporedne (m1) (" << at10->chain_id << "," << at10->resn << ") (" << at11->chain_id << "," 
           << at11->resn << ") (" << at12->chain_id << "," << at12->resn << ")");
#endif

//        float min_dist = 999999.9;
      float max_S = -999999.9;

      vector<Atom*> trojka;

      trojka.clear(); 

      float S = 0.0;

      /* 5. v radiju NEKAJ angstremov od centra trojke iscemo trojko na m2 */ 
      for (map<ChResi, Atom*>::iterator it2 = ca2.begin(); it2 != --(--ca2.end()); it2++) {
        Atom *at20, *at21, *at22;

        map<ChResi, Atom*>::iterator hit2 = it2;
        at20 = hit2->second;

        hit2++; at21 = hit2->second;
        hit2++; at22 = hit2->second;

        /* dobimo najblizjo trojko treh zaporednih aminokislin, vse pripadajo isti verigi na molekuli 2 */
        if ((at20->chain_id == at21->chain_id) && (at21->chain_id == at22->chain_id) && (at20->resi == at21->resi - 1) && (at21->resi == at22->resi - 1)) {

          float d = dist(at11->crd, at21->crd);
          
#ifndef NDEBUG
          dbgmsg("Subgraph::extend_alignment Distance = " << d );
          dbgmsg("Subgraph::extend_alignment Tri zaporedne (m2) (" << at20->chain_id << "," << at20->resn << ") (" << at21->chain_id << "," 
               << at21->resn << ") (" << at22->chain_id << "," << at22->resn << ")");
#endif

          /* le za tiste trojke, ki so si blizu izracunamo S */
          if (d < GLOBAL_DIST) {
            
            S = score->score_blosum(one_letter_code(at10->resn), one_letter_code(at20->resn)) +
              score->score_blosum(one_letter_code(at11->resn), one_letter_code(at21->resn)) +
              score->score_blosum(one_letter_code(at12->resn), one_letter_code(at22->resn));

            if(S > max_S && S > GLOBAL_SVAL) {
#ifndef NDEBUG
              dbgmsg("Subgraph::extend_alignment Distance = " << d );
#endif
              max_S = S;

              trojka.clear();

              trojka.push_back(at20);
              trojka.push_back(at21);
              trojka.push_back(at22);
            }
          }
          
        }
      }

      /* preverimo ali sta trojki na obeh verigah priblizno vzporedni */

      if (!trojka.empty()) {

        ChResi l10(at10->chain_id, at10->resi);
        ChResi l11(at10->chain_id, at10->resi+1);
        ChResi l12(at10->chain_id, at10->resi+2);
        ChResi l13(at10->chain_id, at10->resi+3);
        ChResi l14(at10->chain_id, at10->resi+4);
        ChResi r11(at10->chain_id, at10->resi-1);
        ChResi r12(at10->chain_id, at10->resi-2);

        ChResi l20(trojka[0]->chain_id, trojka[0]->resi);
        ChResi l21(trojka[0]->chain_id, trojka[0]->resi+1);
        ChResi l22(trojka[0]->chain_id, trojka[0]->resi+2);
        ChResi l23(trojka[0]->chain_id, trojka[0]->resi+3);
        ChResi l24(trojka[0]->chain_id, trojka[0]->resi+4);
        ChResi r21(trojka[0]->chain_id, trojka[0]->resi-1);
        ChResi r22(trojka[0]->chain_id, trojka[0]->resi-2);


        Coor vector1, vector2;

        set_zero(vector1);
        set_zero(vector2);

        /* upostevamo (poleg trojke) se dva dodatna C-alpha na vsaki strani (tako dobimo daljsa vektorja) */
        if (ca1.find(l13) != ca1.end() && ca1.find(r11) != ca1.end())
          vector1 = ca1[l12]->crd - ca1[l10]->crd + ca1[l13]->crd - ca1[r11]->crd;

        /* ce manjka l13 ali l23 racunamo vektor z eno ak bolj v levo */
        else if (ca1.find(l13) == ca1.end() && ca1.find(r11) != ca1.end() && ca1.find(r12) != ca1.end())
          vector1 = ca1[l11]->crd - ca1[r11]->crd + ca1[l12]->crd - ca1[r12]->crd;

        /* ce manjka r11 ali r21 se premaknemo za eno v desno */
        else if (ca1.find(r11) == ca1.end() && ca1.find(l13) != ca1.end() && ca1.find(l14) != ca1.end())
          vector1 = ca1[l13]->crd - ca1[l11]->crd + ca1[l14]->crd - ca1[l10]->crd;





        if (ca2.find(l23) != ca2.end() && ca2.find(r21) != ca2.end())
          vector2 = ca2[l22]->crd - ca2[l20]->crd + ca2[l23]->crd - ca2[r21]->crd;

        else if (ca2.find(l23) == ca2.end() && ca2.find(r21) != ca2.end() && ca2.find(r22) != ca2.end())
          vector2 = ca2[l21]->crd - ca2[r21]->crd + ca2[l22]->crd - ca2[r22]->crd;

        else if (ca2.find(r21) == ca2.end() && ca2.find(l23) != ca2.end() && ca2.find(l24) != ca2.end())
          vector2 = ca2[l23]->crd - ca2[l21]->crd + ca2[l24]->crd - ca2[l20]->crd;

        /* ce enega ali obeh vektorjev nismo mogli dolociti, avtomatsko zbrisemo trojko*/
        if (vector1 == null_vect || vector2 == null_vect) {
          trojka.clear();
        }
        else {
          /* v nasprotem primeru pa zbrisemo trojko samo, ce je kot med vektorjema prevelik */
          float ang = angle(vector1, vector2, vector1 % vector2);
          
#ifndef NDEBUG
          dbgmsg("Subgraph::extend_alignment DELETED Angle between [" << at12->resn << at12->resi << at12->chain_id << ","
               << at10->resn << at10->resi << at10->chain_id << "] and ["
               << trojka[2]->resn << trojka[2]->resi << trojka[2]->chain_id << ","
               << trojka[0]->resn << trojka[0]->resi << trojka[0]->chain_id << "] is "
               << ang);
#endif
          
          if (ang > GLOBAL_ANGLE) {
            
            trojka.clear();
          }
        }

      }


      /* 6. ce trojko najdemo, izracunamo score, sicer nazaj na tocko 4 */
      if (!trojka.empty()) {
        /* 7. ce je ta score boljsi od NEKEGA cutoffa, potem gremo naprej, sicer nazaj na tocko 4 */

#ifndef NDEBUG
        dbgmsg("Subgraph::extend_alignment S-val (trojka) = " << max_S );
        //          dbgmsg("Subgraph::extend_alignment E-val (trojka) = " << E );
        dbgmsg("Subgraph::extend_alignment Dobra trojka (m1) " << at10->resn << at10->resi << at10->chain_id << "," 
             << at11->resn << at11->resi << at11->chain_id << "," 
             << at12->resn << at12->resi << at12->chain_id);
        dbgmsg("Subgraph::extend_alignment Dobra trojka (m2) " << trojka[0]->resn << trojka[0]->resi << trojka[0]->chain_id << "," 
             << trojka[1]->resn << trojka[1]->resi << trojka[1]->chain_id << "," 
             << trojka[2]->resn << trojka[2]->resi << trojka[2]->chain_id);
#endif

          /* 8. odpremo nov HSSP in dodamo obe trojki */
//          if (S > GLOBAL_SVAL) HSSP.push_back(make_pair(make_pair(at10->resi, trojka[0]->resi), 3));
        HSSP.push_back(make_pair(TwoChResi(ChResi(at10->chain_id, at10->resi), ChResi(trojka[0]->chain_id, trojka[0]->resi)), 3));
      }
    }
  }

  /* 9. podaljsujemo vse HSSP v levo in desno brez gapov, dokler je njihov score > NEKEGA cutoffa */
  for (vector < pair <TwoChResi, int > >::iterator it = HSSP.begin(); it != HSSP.end(); it++) {

#ifndef NDEBUG
    dbgmsg("Subgraph::extend_alignment Dolzina HSSP = " << it->second );
#endif

    ChResi cr1 = it->first.first;
    ChResi cr2 = it->first.second;

    char ch1 = cr1.first;
    char ch2 = cr2.first;
    int i1 = cr1.second;
    int i2 = cr2.second;


    /* najprej v desno */
    float S_prvi = score->score_blosum(one_letter_code(ca1[ChResi(ch1, i1 + 0)]->resn), one_letter_code(ca2[ChResi(ch2, i2 + 0)]->resn)) +
                   score->score_blosum(one_letter_code(ca1[ChResi(ch1, i1 + 1)]->resn), one_letter_code(ca2[ChResi(ch2, i2 + 1)]->resn)) +
                   score->score_blosum(one_letter_code(ca1[ChResi(ch1, i1 + 2)]->resn), one_letter_code(ca2[ChResi(ch2, i2 + 2)]->resn));


    float E_prvi = K * m1->sequence.size() * m2->sequence.size() * exp( (-1) * LAMBDA * S_prvi );

    float S = S_prvi;
    float E = E_prvi;

    float E_prej = E_prvi;

    int offs = 3;

    /* sele od offset je 3 upostevamo tudi E-VAL (tako sigurno dobimo trojko v HSSP) */
    while (ca1.find(ChResi(ch1, i1 + offs)) != ca1.end() && ca2.find(ChResi(ch2, i2 + offs)) != ca2.end()) {

      S += score->score_blosum(one_letter_code(ca1[ChResi(ch1, i1 + offs)]->resn), one_letter_code(ca2[ChResi(ch2, i2 + offs)]->resn));
      E_prej = E;

      E = K * m1->sequence.size() * m2->sequence.size() * exp( (-1) * LAMBDA * S );

      if (E < E_prej) offs++; else break;

    }
    /* nova dolzina HSSP */ 
    it->second = offs ;

#ifndef NDEBUG
    dbgmsg("Subgraph::extend_alignment Dolzina HSSP = " << it->second );
#endif
    /* nato se v levo */
    S = S_prvi;
    E = E_prvi;
    E_prej = E_prvi;

    offs = -1;

    while (ca1.find(ChResi(ch1, i1 + offs)) != ca1.end() && ca2.find(ChResi(ch2, i2 + offs)) != ca2.end()) {

      S += score->score_blosum(one_letter_code(ca1[ChResi(ch1, i1 + offs)]->resn), one_letter_code(ca2[ChResi(ch2, i2 + offs)]->resn));
      E_prej = E;

      E = K * m1->sequence.size() * m2->sequence.size() * exp( (-1) * LAMBDA * S );
      if (E < E_prej) offs--; else break;
      
    }
    /* nova dolzina in novi zacetki HSSP */
    it->second -= (offs + 1);
    it->first.first.second += (offs + 1);
    it->first.second.second += (offs + 1);

#ifndef NDEBUG
    dbgmsg("Subgraph::extend_alignment Dolzina HSSP = " << it->second );
#endif
  }
  /*10. sortiramo od najvecjega do najmanjsega HSSP-ja */
  stable_sort(HSSP.begin(), HSSP.end(), by_hssp_length);

  /*11. izberemo HSSP-je, ki se ne prekrivajo */
  for (vector < pair < TwoChResi, int > >::iterator it = HSSP.begin(); it != HSSP.end(); it++) {

    bool prekrivata = false;
#ifndef NDEBUG
    dbgmsg("Subgraph::extend_alignment Dolzina HSSP (po sortiranju) = " << it->second );
#endif

    for (vector < pair < TwoChResi, int > >::iterator it2 = HSSP.begin(); it2 != it; it2++) {
      /* dolzina HSSP != -1 pomeni, da je bil dodan v kliko */
      if (it->second != -1) {

        int first_i1 = it->first.first.second;
        int first_i2 = it->first.second.second;
        int last_i1 = first_i1 + it->second;
        int last_i2 = first_i2 + it->second;
        
        int first_j1 = it2->first.first.second;
        int first_j2 = it2->first.second.second;
        int last_j1 = first_j1 + it2->second;
        int last_j2 = first_j2 + it2->second;
        
        if ((first_i1 >= first_j1 && first_i1 < last_j1) || (last_i1 >= first_j1 && last_i1 < last_j1) ||
            (first_i2 >= first_j2 && first_i2 < last_j2) || (last_i2 >= first_j2 && last_i2 < last_j2))
          prekrivata = true;
      }
    }
    /*12. ce se HSSP, ki ga predstavlja it, ne prekriva z nobenim ze dodanim, potem ga dodamo v kliko */
    if (!prekrivata) {
      int hssp_length = it->second;

      ChResi cr1 = it->first.first;
      ChResi cr2 = it->first.second;


#ifndef NDEBUG
      dbgmsg("Subgraph::extend_alignment V kliko dodajamo HSSP dolzine = " << it->second );
#endif      
      for (int j = 0; j < hssp_length; j++) {
        
        /* Row r mora biti unique za vsak par prileganih aminokislin */
        map<TwoChResi, Row*>::iterator it2 = row_map.find(TwoChResi(ChResi(cr1.first, cr1.second + j), ChResi(cr2.first, cr2.second + j)));
        if (it2 == row_map.end() ) {
          /* dodaj row v product graf, da bos na koncu sprostil pomnilnik */
          Row *r = new Row(); 
          r->next = p->row;
          p->row = r;

          /* dodaj nova deskriptorja med desc, iz istega razloga */
          Descriptor *d1 = new Descriptor(-1); // -1 oznacuje slabo strukturno ohranjenost
          d1->next = desc1;
          desc1 = d1;
          Descriptor *d2 = new Descriptor(-1);
          d2->next = desc2;
          desc2 = d2;

          d1->atom = ca1[ChResi(cr1.first, cr1.second + j)];
          d2->atom = ca2[ChResi(cr2.first, cr2.second + j)];

#ifndef NDEBUG
          dbgmsg("("<< one_letter_code(d1->atom->resn) << d1->atom->resi << d1->atom->chain_id << ","
               << one_letter_code(d2->atom->resn) << d2->atom->resi << d2->atom->chain_id << ")");
#endif
          r->desc1 = d1;
          r->desc2 = d2;

          clus[clNum]->vert.push_back(new Vert(0, 0, 0, r ) );

          row_map[TwoChResi(ChResi(cr1.first, cr1.second + j), ChResi(cr2.first, cr2.second + j))] = r;
        }
        else {
          clus[clNum]->vert.push_back(new Vert(0, 0, 0, it2->second ) );
        }
      }
#ifndef NDEBUG
      dbgmsg(endl);
#endif
    }
    else {
      /* oznacimo, da tega HSSP-ja nismo dodali v kliko */
      it->second = -1;
    }
  }


  /* 13. nazaj rotiramo m1 na prvotne koordinate, tako kot narekuje ROTA in TRANS v clus[clNum]  */
  for (Atom *atm = m1->atom; atm != NULL; atm = atm->next) {
    
    /* upostevamo le backbone C-alpha (CA) atome obeh proteinov */
    if (!strcmp(atm->tag, "CA") ) {

#ifndef NDEBUG
      dbgmsg("Subgraph::extend_alignment Koordinate pred rotacijo " << atm->crd.x << "," << atm->crd.y << "," << atm->crd.z );
#endif
      Coor crot = kabsch->rotate_vector_inv(atm->crd, clus[clNum]->U, clus[clNum]->t);
        
      atm->crd = crot;
      
#ifndef NDEBUG
      dbgmsg("Subgraph::extend_alignment Koordinate po rotaciji " << atm->crd.x << "," << atm->crd.y << "," << atm->crd.z );
#endif

//      ca1[atm->resi] = atm;

    }
  }
  

  delete kabsch;

}



void Subgraph::output_qmax(Output *out, Score *score) {
  for (unsigned int i = 0; i < qmax.size(); i++) {
      dbgmsg("START> CLIQUE " << i );
      out->out_score(qmax[i]);

//~ #ifndef NDEBUG
      //~ out->out_desc(qmax[i]);
      //~ out->out_resi(qmax[i]);
      //~ out->output_align_protein(mol1, mol2, qmax[i]);
//~ #endif
#ifndef VERB
      out->out_desc_xquiet(qmax[i]);
      out->out_resi_quiet(qmax[i], score);
#endif
      out->out_rota(qmax[i]->U);
      out->out_trans(qmax[i]->t);
      dbgmsg("END> CLIQUE" );
  }
}

void Subgraph::output_clus(Output *out, Score *score) {
  for (unsigned int i = 0; i < clus.size(); i++) {
      dbgmsg("START> CLUSTER " << i );
      out->out_score(clus[i]);

//~ #ifndef NDEBUG
      //~ out->out_desc(clus[i]);
      //~ out->out_resi(clus[i]);
      //~ out->output_align_protein(mol1, mol2, clus[i]);
//~ #endif
#ifndef VERB
      out->out_desc_xquiet(clus[i]);
      out->out_resi_quiet(clus[i], score);
#endif
      out->out_rota(clus[i]->U);
      out->out_trans(clus[i]->t);
      dbgmsg("END> CLUSTER" );
  }
}


//void Subgraph::delete_bad_scoring() {
//  /* 
//     Delete all items from qmax array that have unfeasible surface complementarities or
//     bad blosum scores. 
//     NEW (AUG/19/2009) : izbrisi tudi vse s slabim rmsd-jem med C-alpha atomi :)
//     NOVO (APR/7/2010) : ce je calpha_rmsd = -1 , potem je prislo do napake pri superimponiranju (superimpose_calpha vrne -1)
//                         to naredimo zato, ker pri racunanju cluster_score, ne sme biti negativne vrednosti za calpha_rmsd
//  */
//  dbgmsg("Deleting bad scoring cliques..." );
//  int num = 0;
//  for (unsigned int i = 0; i < qmax.size(); i++) {
//    if (z_score(qmax[i]->cluster_score) < Z_SCORE ||
//        qmax[i]->calpha_rmsd < 0) {
//
//      delete qmax[i];
//      qmax[i] = qmax.back();
//      
//      qmax.pop_back();
//      i--; 
//      num++;
//    }
//  }
//  dbgmsg("Number of bad scoring cliques deleted = " << num );
//}

void Subgraph::delete_bad_scoring() {
  /* 
     Delete all items from qmax array that have unfeasible surface complementarities or
     bad blosum scores. 
     NEW (AUG/19/2009) : izbrisi tudi vse s slabim rmsd-jem med C-alpha atomi :)
     NOVO (APR/7/2010) : ce je calpha_rmsd = -1 , potem je prislo do napake pri superimponiranju (superimpose_calpha vrne -1)
                         to naredimo zato, ker pri racunanju cluster_score, ne sme biti negativne vrednosti za calpha_rmsd
  */
  dbgmsg("Deleting bad scoring cliques..." );
  int num = 0;
  for (unsigned int i = 0; i < qmax.size(); i++) {
    if (qmax[i]->surf_vector_angle > SURF_VECTOR_ANGLE ||
        qmax[i]->calpha_rmsd > CALPHA_RMSD || 
        qmax[i]->calpha_rmsd < 0 || 
//        qmax[i]->vert.size() < (unsigned) SCONS ||
        qmax[i]->blosum_score > BLOSUM_SCORE) {

      delete qmax[i];
      qmax[i] = qmax.back();
      
      qmax.pop_back();
      i--; 
      num++;
    }
  }
  dbgmsg("Number of bad scoring cliques deleted = " << num );
}

void Subgraph::sort_qmax(bool (*callback)(Item *i, Item *j)) {
  stable_sort(qmax.begin(), qmax.end(), callback);
}

void Subgraph::sort_clus(bool (*callback)(Item *i, Item *j)) {
  stable_sort(clus.begin(), clus.end(), callback);
}

void Subgraph::sort_item(bool (*callback)(Item *i, Item *j)) {
  stable_sort(item.begin(), item.end(), callback);
}



Subgraph::~Subgraph() {
  dbgmsg("Deleting " << item.size() << " subgraphs ... " );
  dbgmsg("Deleting " << qmax.size() << " cliques ... " );
  dbgmsg("Deleting " << clus.size() << " clusters ... " );
  for (unsigned int i = 0; i < item.size(); i++) 
    delete item[i];
  for (unsigned int i = 0; i < qmax.size(); i++) 
    delete qmax[i];
  for (unsigned int i = 0; i < clus.size(); i++) 
    delete clus[i];

  item.clear();
  qmax.clear();
  clus.clear();
  dbgmsg("Done" );
}

Subgraph::Subgraph() {
  /*
    Another constructor creating just an empty subgraph object - don't intialize with product graph.
    We need this for --multi option, where subgraphs are being read from the mysql database.
    NOTE: we do not initialize Item *qmax nor Item *item here but in the read_graph function of mysql.cc
  */
//  item = NULL; qmax = NULL; clus = NULL; size = 0;


}

Subgraph::Subgraph(Product *p) {
  /*
    We read all subgraphs at once into array subgraph. In item[r->sgraph].vert[ i ] are 
    all vertices of this sgraph and in item[r->sgraph].edge[ j ] are all connected pairs of vertices.
    The vertex vert[ 0 ] is connected to all the vertices in the subgraph, that's why it has 
    to be dealt with separately.
    NOTE: we do not initialize Item *clus here but in the join function of cluster.cc
  */
  dbgmsg("Read subgraphs from product graph (enhanced)... " );
  Row *r;
  Column *c, *c2;
  double dist1, dist2;
  int i, j, k;
  Column *temp_c[2000];

//  item = NULL; qmax = NULL; clus = NULL; size = 0;

  item.reserve(p->num_subgraph);
  qmax.reserve(p->num_subgraph);

//  size = p->num_subgraph;
//  item = new Item[size];
//  qmax = new Item[size];

//  item = (Item*) calloc(size, sizeof(Item));
//  qmax = (Item*) calloc(size, sizeof(Item));

#ifndef NDEBUG
  dbgmsg("Number of subgraphs = " << p->num_subgraph );
#endif
#ifndef NDEBUG
  clock_t start = clock();
#endif
  r = p->row;
  while (r->next != NULL) {

//    if (r->score > MNSP_HIGH) { 
#ifndef NDEBUG
      dbgmsg("SUBGRAPH NUMBER = " << r->sgraph );
#endif
//      qmax[r->sgraph].U = gsl_matrix_alloc(3, 3);
//      qmax[r->sgraph].t = gsl_vector_alloc(3);
      c = r->column; i = 0;
      while (c->next) { i++; c = c->next; } // c->num begins at 1
#ifndef NDEBUG
      dbgmsg("SUBGRAPH SIZE = " << i+1 );
#endif
//      item[r->sgraph].size = i+1;
      item.push_back(new Item() );
      qmax.push_back(new Item() );
      
      item.back()->vert.reserve( i + 1 );
      qmax.back()->vert.reserve( i + 100 );  // naj bo qmax.vert vecji (zaradi dodatnih tock) !!!

//      item[r->sgraph].vert.reserve( i + 1 );
//      qmax[r->sgraph].vert.reserve( i + 100 );  // naj bo qmax.vert vecji (zaradi dodatnih tock) !!!

//      item[r->sgraph].vert = (Vert*) calloc(i+1, sizeof(Vert));
//      qmax[r->sgraph].vert = (Vert*) calloc(i+100, sizeof(Vert));  // naj bo qmax.vert vecji (zaradi dodatnih tock) !!!

      item.back()->edge.reserve( (i + 1) * i / 2 ); // to minimize reallocations, we expect the max number of edges 

//      item[r->sgraph].edge.reserve( (i + 1) * i / 2 ); // to minimize reallocations, we expect the max number of edges 

      item.back()->vert.push_back(new Vert(0, 0, i, r ) );  // leading vertex 

//      item[r->sgraph].vert.push_back(new Vert(0, 0, i, r ) );  // leading vertex 

//      item[r->sgraph].vert[ 0 ][ 0 ].vertex = 0;
//      item[r->sgraph].vert[ 0 ][ 0 ].row = r;
      
//      item[r->sgraph].vert[0].vertex = 0;
//      item[r->sgraph].vert[0].row = r;
//      item[r->sgraph].vert[0].item.size = i;
//      item[r->sgraph].vert[0].item.vert = (Vert*) calloc(i, sizeof(Vert));
     
      c = r->column;
      while (c->next) {

        item.back()->edge.push_back(make_pair(0, c->num ) );  // insert edges for the leading vertex

//        item[r->sgraph].edge.push_back(make_pair(0, c->num ) );  // insert edges for the leading vertex

//        item[r->sgraph].vert[ 0 ].push_back(new Vertex(c->num ) );  // it's connected vert, no need to specify row

//        item[r->sgraph].vert[ 0 ][ c->num ].vertex = c->num;
//        item[r->sgraph].vert[0].item.vert[c->num - 1].vertex = c->num;
        c = c->next;
      }

      c = r->column;
      while (c->next) {


//        item[r->sgraph].vert[ c->num ][ 0 ].vertex = c->num;
//        item[r->sgraph].vert[ c->num ][ 0 ].row = c->item;

//        item[r->sgraph].vert[c->num].vertex = c->num;
//        item[r->sgraph].vert[c->num].row = c->item;

        j = 0;
        c2 = r->column;
        while (c2->next) {
          if (c != c2) {
            dist1 = dist_fast(c->item->desc1->crd, c2->item->desc1->crd );
            dist2 = dist_fast(c->item->desc2->crd, c2->item->desc2->crd );
            if (pow(dist1 + dist2 - RESOLUTION*RESOLUTION, 2 ) < 4*dist1*dist2 )
              if ( (sqrt(dist1) + sqrt(dist2) ) / 2 < CUTOFF_SECOND )
                //              if (sqrt(dist1) < CUTOFF_SECOND && sqrt(dist2) < CUTOFF_SECOND)
                //              if (sqrt(dist1) < CUTOFF_SECOND)
                temp_c[ j++ ] = c2;
          }       
          c2 = c2->next;
        }

        item.back()->vert.push_back(new Vert(c->num, 0, j + 1, c->item ) );

//        item[r->sgraph].vert.push_back(new Vert(c->num, 0, j + 1, c->item ) );


//        item[r->sgraph].vert[ c->num ].reserve( j + 2 );
//        item[r->sgraph].vert[ c->num ][ 1 ].vertex = 0;

//        item[r->sgraph].vert[c->num].item.size = j+1;
//        item[r->sgraph].vert[c->num].item.vert = (Vert*) calloc(j+1, sizeof(Vert));
//        item[r->sgraph].vert[c->num].item.vert[0].vertex = 0;

        item.back()->edge.push_back(make_pair(c->num, 0 ) );

//        item[r->sgraph].edge.push_back(make_pair(c->num, 0 ) );

        k = j+1;
        j = 1;
        while (j < k) {
          item.back()->edge.push_back(make_pair(c->num, temp_c[ j - 1 ]->num ) );

//          item[r->sgraph].edge.push_back(make_pair(c->num, temp_c[ j - 1 ]->num ) );

//          item[r->sgraph].vert[ c->num ][ j + 1 ].vertex = temp_c[ j - 1 ]->num;

//          item[r->sgraph].vert[ c->num ].item.vert[j].vertex = temp_c[j-1]->num;
          j++;
          //          dbgmsg("j = " << j );
        }            
        c = c->next;
      }
#ifndef NDEBUG
      dbgmsg("VERB> SUBGRAPH " << r->sgraph );
      dbgmsg("VERB> Vertices " );
      for (unsigned int i = 0; i < item.back()->vert.size(); i++)
        dbgmsg("VERB> v=" << item.back()->vert[i]->vertex << 
          " se=" << item.back()->vert[i]->size_edge);
      dbgmsg("VERB> Edges " );
      for (unsigned int i = 0; i < item.back()->edge.size(); i++)
        dbgmsg("VERB> " << item.back()->edge[i].first << 
          " " << item.back()->edge[i].second);
#endif
//    }
    r = r->next;
  }
  dbgmsg("Subgraph::Subgraph Time = " << (double)(clock() - start) / CLOCKS_PER_SEC );
//  exit(1);
}




