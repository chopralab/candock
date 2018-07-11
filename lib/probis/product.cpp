#include "product.h"
#include "atom.h"
#include "debug.hpp"
#include "desc.h"
#include "eelement.h"
#include "item.h"
#include "score.h"

bool by_mnsp(Descriptor* di, Descriptor* dj) {
    return (di->s->mnsp < dj->s->mnsp);
}  // callback for sorting vertices by degree
bool by_desc_num(Descriptor* di, Descriptor* dj) {
    return (di->num < dj->num);
}  // callback for sorting vertices by degree

Product::~Product() {
    dbgmsg("Deleting product graph ...");
    Row* r;
    Column* c;
    while (row != NULL) {
        r = row;
        while (row->column != NULL) {
            c = row->column;
            row->column = row->column->next;

            delete c->item;  // ker so verteksi product grafa (ki so tipa Row*),
                             // pointerji, ki so posebni za vsak podgraf
            delete c;
        }
        row = row->next;
        delete r;
    }

    /* pobrisemo tudi item iz filtrat */
    for (vector<pair<pair<Oglisca, Oglisca>, Item*> >::iterator it =
             filtrat.begin();
         it != filtrat.end(); it++) {
        delete it->second->vert[0]->row;
        delete it->second->vert[1]->row;
        delete it->second->vert[2]->row;
        delete it->second->vert[3]->row;
        delete it->second;
    }

    filtrat.clear();
    thedron1.clear();
    thedron2.clear();
}

void Product::tetrahedron(Score* s, Descriptor* desc1, Descriptor* desc2,
                          int psurf) {
    /*
      Poiscemo podobne tetraedre v obeh proteinih, pri cemer so oglisca
      deskriptorji z enakimi oznakami.
      Note: psurf = 0 when total proteins search, psurf = 1 when protein-protein
      interaction search
    */

    dbgmsg("Searching for tetrahedrons ...");
#ifndef NDEBUG
    clock_t start = clock();
#endif
    vector<float> d;
    d.reserve(6);

    vector<Descriptor*> de;
    de.reserve(4);

    /* bitne maske */
    int ace = 0X1000;
    int don = 0X200;
    int ado = 0X40;
    int ali = 0X08;
    int pic = 0X01;

    thedron1.clear();
    thedron2.clear();

    /* najdemo vse tetraedre v prvem proteinu */
    Descriptor* d1 = desc1;

    while (d1 != NULL) {
        if (d1->psurf == psurf &&
            !d1->bb) {  // pri PPI search ima psurf 1 binding site (interface),
                        // s katerimi searchas (funkcija trim)

            //      dbgmsg("deskriptorji = " << d1->atom->resn << d1->atom->resi
            //      << " mnsp = " << d1->s->mnsp);
            /* po vseh sosedih deskriptorja */
            for (Element* sosed1 = d1->neighb; sosed1 != NULL;
                 sosed1 = sosed1->next) {
                /* ce obiscemo kaksen predhoden (ze obravnavan) deskriptor, ni
                 * treba nadaljevati gradnje tetraedra */
                Descriptor* d2 = (Descriptor*)sosed1->item;
                if (d2->num < d1->num &&
                    !d2->bb) {  // num tecejo od vecje k manjsi

                    for (Element* sosed2 = sosed1->next; sosed2 != NULL;
                         sosed2 = sosed2->next) {
                        Descriptor* d3 = (Descriptor*)sosed2->item;
                        if (d3->num < d1->num && !d3->bb) {
                            for (Element* sosed3 = sosed2->next; sosed3 != NULL;
                                 sosed3 = sosed3->next) {
                                Descriptor* d4 = (Descriptor*)sosed3->item;

                                /* deskriptorji morajo izhajati vsaj iz treh
                                 * razlicnih aminokislin */
                                /* prave kombinacije so 123x 12x3 1x23 x123 (x
                                 * je lahko tudi enaka ostalim trem) */

                                // z multi chain ne moremo vec uporabiti
                                // atom->resi za razlikovanje med aminokislinami
                                if (d4->num < d1->num && !d4->bb &&
                                    (((d1->atom->acid_start !=
                                       d2->atom->acid_start) &&
                                      ((d1->atom->acid_start !=
                                            d3->atom->acid_start &&
                                        d2->atom->acid_start !=
                                            d3->atom->acid_start) ||
                                       (d1->atom->acid_start !=
                                            d4->atom->acid_start &&
                                        d2->atom->acid_start !=
                                            d4->atom->acid_start))) ||
                                     ((d3->atom->acid_start !=
                                       d4->atom->acid_start) &&
                                      ((d1->atom->acid_start !=
                                            d3->atom->acid_start &&
                                        d1->atom->acid_start !=
                                            d4->atom->acid_start) ||
                                       (d2->atom->acid_start !=
                                            d3->atom->acid_start &&
                                        d2->atom->acid_start !=
                                            d4->atom->acid_start))))) {
                                    /* v kodo tudi dolzine najvecjih treh
                                     * stranic, sortirane od najvecje do
                                     * najmanjse */

                                    d[0] = dist(d1->crd, d2->crd);
                                    d[1] = dist(d1->crd, d3->crd);
                                    d[2] = dist(d1->crd, d4->crd);
                                    d[3] = dist(d2->crd, d3->crd);
                                    d[4] = dist(d2->crd, d4->crd);
                                    d[5] = dist(d3->crd, d4->crd);

                                    //                  sort(d.begin(),d.end());
                                    sort(d.begin(), d.begin() + 6);

                                    /* ce je najvecja stranica manjsa od
                                     * CUTOFF_FIRST (=6A) */
                                    if (d[5] < CUTOFF_FIRST) {
                                        long int koda = 0;

                                        //                  dbgmsg("tu sem " <<
                                        //                  d[5] << " " << d[4]
                                        //                  << " " << d[3] << "
                                        //                  " << d[2] << " " <<
                                        //                  d[1] << " " << d[0]
                                        //                  << " " );
                                        //                  for
                                        //                  (vector<float>::iterator
                                        //                  dit = d.begin(); dit
                                        //                  != d.end(); dit++)
                                        //                  dbgmsg(thedron1.size()
                                        //                  << " " << *dit);

                                        //                    koda += (long int)
                                        //                    floor((d[5] +
                                        //                    0.5)/2) *
                                        //                    0X200000;
                                        //                    koda += (long int)
                                        //                    floor((d[4] +
                                        //                    0.5)/2) *
                                        //                    0X1000000;
                                        //                    koda += (long int)
                                        //                    floor((d[3] +
                                        //                    0.5)/2) *
                                        //                    0X8000000;

                                        //                  koda += floor(d[2] +
                                        //                  0.5) * 0X40000000;
                                        //                  koda += floor(d[1] +
                                        //                  0.5) * 0X200000000;
                                        //                  koda += floor(d[0] +
                                        //                  0.5) * 0X1000000000;

                                        /*  koda iz mnsp deskriptorjev - to bo
                                         * key za hash tabelo ace|don|ado|al|pi
                                         * vsak ima po tri bite */
                                        switch (d1->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }
                                        switch (d2->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }
                                        switch (d3->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }
                                        switch (d4->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }

                                        /* v kodo dodamo tudi volumen tetraedra
                                         */

                                        float V =
                                            fabs(((d1->crd - d4->crd) *
                                                  ((d2->crd - d4->crd) %
                                                   (d3->crd - d4->crd)))) /
                                            6;

                                        //                  dbgmsg(" V = " << V
                                        //                  << " round(V) = " <<
                                        //                  floor(V + 0.5));
                                        //                    koda += (long int)
                                        //                    floor(V/2 + 0.5) *
                                        //                    0X8000;

                                        Oglisca og;
                                        //                  og.d1 = d1;og.d2 =
                                        //                  d2;og.d3 = d3;og.d4
                                        //                  = d4;

                                        de[0] = d1;
                                        de[1] = d2;
                                        de[2] = d3;
                                        de[3] = d4;

                                        //                    sort(de.begin(),de.end(),
                                        //                    by_mnsp);
                                        //                    sort(de.begin(),de.begin()
                                        //                    + 4, by_mnsp);

                                        //                    dbgmsg("tu sem "
                                        //                    << de[0]->s->mnsp
                                        //                    << " " <<
                                        //                    de[1]->s->mnsp <<
                                        //                    " " <<
                                        //                    de[2]->s->mnsp <<
                                        //                    " " <<
                                        //                    de[3]->s->mnsp);

                                        /* sortiramo jih po stevilki
                                         * deskriptorja desc->num */
                                        sort(de.begin(), de.begin() + 4,
                                             by_desc_num);

                                        //                    cout <<dec<< "1 "
                                        //                         <<
                                        //                         de[0]->atom->resn
                                        //                         <<de[0]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[1]->atom->resn
                                        //                         <<de[1]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[2]->atom->resn
                                        //                         <<de[2]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[3]->atom->resn
                                        //                         <<de[3]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[0]->s->mnsp
                                        //                         << " "
                                        //                         <<
                                        //                         de[1]->s->mnsp
                                        //                         << " "
                                        //                         <<
                                        //                         de[2]->s->mnsp
                                        //                         << " "
                                        //                         <<
                                        //                         de[3]->s->mnsp
                                        //                         << " "
                                        //                         << " V = " <<
                                        //                         V << " d[5] =
                                        //                         " << d[5] <<
                                        //                         " d[4] = " <<
                                        //                         d[4] << "
                                        //                         d[3] = " <<
                                        //                         d[3]
                                        //                         << " d[2] = "
                                        //                         << d[2] << "
                                        //                         d[1] = " <<
                                        //                         d[1] << "
                                        //                         d[0] = " <<
                                        //                         d[0]<<endl;
                                        //                    cout << " koda = "
                                        //                    <<oct<< koda
                                        //                         << endl;

                                        /* najdemo vse permutacije deskriptorjev
                                         * ( 1a 1b 4 8 -> 1b 1a 4 8 ) samo
                                         * deskriptorje z enakimi mnsp
                                         * permutiramo */
                                        /* to je treba storiti le za thedron1!
                                         * za thedron2 ni treba!! */
                                        do {
                                            /* zapisemo v thedron1 le
                                             * permutacije, ki imajo mnsp-je v
                                             * narascajocem vrstnem redu (mora
                                             * biti enako kot by_mnsp funkcija)
                                             */
                                            if (de[0]->s->mnsp <=
                                                    de[1]->s->mnsp &&
                                                de[1]->s->mnsp <=
                                                    de[2]->s->mnsp &&
                                                de[2]->s->mnsp <=
                                                    de[3]->s->mnsp) {
                                                og.d1 = de[0];
                                                og.d2 = de[1];
                                                og.d3 = de[2];
                                                og.d4 = de[3];

                                                og.a = d[5];
                                                og.b = d[4];
                                                og.c = d[3];

                                                og.V = V;

                                                /* dodamo tetraeder v hash
                                                 * tabelo (STL multimap) */
                                                thedron1.insert(
                                                    pair<long int, Oglisca>(
                                                        koda, og));
#ifndef NDEBUG
                                                dbgmsg(de[0]->s->mnsp
                                                       << " " << de[1]->s->mnsp
                                                       << " " << de[2]->s->mnsp
                                                       << " "
                                                       << de[3]->s->mnsp);
                                                dbgmsg(de[0]->num
                                                       << " " << de[1]->num
                                                       << " " << de[2]->num
                                                       << " " << de[3]->num);
#endif
                                            }

                                        } while (next_permutation(
                                            de.begin(), de.begin() + 4,
                                            by_desc_num));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        d1 = d1->next;
    }

    /* najdemo vse tetraedre v drugem proteinu */
    d1 = desc2;

    while (d1 != NULL) {
        if (d1->psurf == 0 && !d1->bb) {  // drugi protein ima psurf vedno 0
                                          // (ker vedno upostevas celega)

            /* po vseh sosedih deskriptorja */
            for (Element* sosed1 = d1->neighb; sosed1 != NULL;
                 sosed1 = sosed1->next) {
                /* ce obiscemo kaksen predhoden (ze obravnavan) deskriptor, ni
                 * treba nadaljevati gradnje tetraedra */
                Descriptor* d2 = (Descriptor*)sosed1->item;
                if (d2->num < d1->num &&
                    !d2->bb) {  // num tecejo od vecje k manjsi

                    for (Element* sosed2 = sosed1->next; sosed2 != NULL;
                         sosed2 = sosed2->next) {
                        Descriptor* d3 = (Descriptor*)sosed2->item;
                        if (d3->num < d1->num && !d3->bb) {
                            for (Element* sosed3 = sosed2->next; sosed3 != NULL;
                                 sosed3 = sosed3->next) {
                                Descriptor* d4 = (Descriptor*)sosed3->item;
                                /* deskriptorji morajo izhajati vsaj iz treh
                                 * razlicnih aminokislin */
                                /* prave kombinacije so 123x 12x3 1x23 x123 (x
                                 * je lahko tudi enaka ostalim trem) */

                                // z multi chain ne moremo vec uporabiti
                                // atom->resi za razlikovanje med aminokislinami
                                if (d4->num < d1->num && !d4->bb &&
                                    (((d1->atom->acid_start !=
                                       d2->atom->acid_start) &&
                                      ((d1->atom->acid_start !=
                                            d3->atom->acid_start &&
                                        d2->atom->acid_start !=
                                            d3->atom->acid_start) ||
                                       (d1->atom->acid_start !=
                                            d4->atom->acid_start &&
                                        d2->atom->acid_start !=
                                            d4->atom->acid_start))) ||
                                     ((d3->atom->acid_start !=
                                       d4->atom->acid_start) &&
                                      ((d1->atom->acid_start !=
                                            d3->atom->acid_start &&
                                        d1->atom->acid_start !=
                                            d4->atom->acid_start) ||
                                       (d2->atom->acid_start !=
                                            d3->atom->acid_start &&
                                        d2->atom->acid_start !=
                                            d4->atom->acid_start))))) {
                                    /* v kodo tudi dolzine najvecjih treh
                                     * stranic, sortirane od najvecje do
                                     * najmanjse */
                                    d[0] = dist(d1->crd, d2->crd);
                                    d[1] = dist(d1->crd, d3->crd);
                                    d[2] = dist(d1->crd, d4->crd);
                                    d[3] = dist(d2->crd, d3->crd);
                                    d[4] = dist(d2->crd, d4->crd);
                                    d[5] = dist(d3->crd, d4->crd);

                                    //                  sort(d.begin(),d.end());
                                    sort(d.begin(), d.begin() + 6);

                                    if (d[5] < CUTOFF_FIRST) {
                                        long int koda = 0;

                                        //                    koda += (long int)
                                        //                    floor(d[5]/2 +
                                        //                    0.5) * 0X200000;
                                        //                    koda += (long int)
                                        //                    floor(d[4]/2 +
                                        //                    0.5) * 0X1000000;
                                        //                    koda += (long int)
                                        //                    floor(d[3]/2 +
                                        //                    0.5) * 0X8000000;

                                        //                  koda += floor(d[2] +
                                        //                  0.5) * 0X40000000;
                                        //                  koda += floor(d[1] +
                                        //                  0.5) * 0X200000000;
                                        //                  koda += floor(d[0] +
                                        //                  0.5) * 0X1000000000;

                                        /*  koda iz mnsp deskriptorjev - to bo
                                         * key za hash tabelo ace|don|ado|al|pi
                                         * vsak ima po tri bite */
                                        switch (d1->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }
                                        switch (d2->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }
                                        switch (d3->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }
                                        switch (d4->s->mnsp) {
                                            case 0X01:
                                                koda += ali;
                                                break;
                                            case 0X02:
                                                koda += pic;
                                                break;
                                            case 0X04:
                                                koda += ace;
                                                break;
                                            case 0X08:
                                                koda += don;
                                                break;
                                            case 0X0C:
                                                koda += ado;
                                                break;
                                        }

                                        /* v kodo dodamo tudi volumen tetraedra
                                         */

                                        float V =
                                            fabs(((d1->crd - d4->crd) *
                                                  ((d2->crd - d4->crd) %
                                                   (d3->crd - d4->crd)))) /
                                            6;

                                        //                  dbgmsg(" V = " << V
                                        //                  << " round(V) = " <<
                                        //                  floor(V + 0.5));
                                        //                    koda += (long int)
                                        //                    floor(V/2 + 0.5) *
                                        //                    0X8000;

                                        Oglisca og;
                                        //                  og.d1 = d1;og.d2 =
                                        //                  d2;og.d3 = d3;og.d4
                                        //                  = d4;

                                        de[0] = d1;
                                        de[1] = d2;
                                        de[2] = d3;
                                        de[3] = d4;

                                        sort(de.begin(), de.begin() + 4,
                                             by_mnsp);

                                        //                    cout <<dec<< "2 "
                                        //                         <<
                                        //                         de[0]->atom->resn
                                        //                         <<de[0]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[1]->atom->resn
                                        //                         <<de[1]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[2]->atom->resn
                                        //                         <<de[2]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[3]->atom->resn
                                        //                         <<de[3]->atom->resi
                                        //                         << " "
                                        //                         <<
                                        //                         de[0]->s->mnsp
                                        //                         << " "
                                        //                         <<
                                        //                         de[1]->s->mnsp
                                        //                         << " "
                                        //                         <<
                                        //                         de[2]->s->mnsp
                                        //                         << " "
                                        //                         <<
                                        //                         de[3]->s->mnsp
                                        //                         << " "
                                        //                         << " V = " <<
                                        //                         V << " d[5] =
                                        //                         " << d[5] <<
                                        //                         " d[4] = " <<
                                        //                         d[4] << "
                                        //                         d[3] = " <<
                                        //                         d[3]
                                        //                         << " d[2] = "
                                        //                         << d[2] << "
                                        //                         d[1] = " <<
                                        //                         d[1] << "
                                        //                         d[0] = " <<
                                        //                         d[0]<<endl;
                                        //                    cout << " koda = "
                                        //                    <<oct<< koda
                                        //                         << endl;

                                        og.d1 = de[0];
                                        og.d2 = de[1];
                                        og.d3 = de[2];
                                        og.d4 = de[3];

                                        og.a = d[5];
                                        og.b = d[4];
                                        og.c = d[3];

                                        og.V = V;

                                        /* dodamo tetraeder v hash tabelo (STL
                                         * multimap) */
                                        thedron2.insert(
                                            pair<long int, Oglisca>(koda, og));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        d1 = d1->next;
    }

    dbgmsg("Product::tetrahedron Vmesni cas = " << (double)(clock() - start) /
                                                       CLOCKS_PER_SEC);

    Kabsch* k = new Kabsch();

    int i = 0;

    filtrat.clear();

    /* gremo po vseh tetraedrih prvega proteina in najdemo korespondence v hash
     * tabeli tetraedrov drugega proteina */
    for (multimap<long int, Oglisca>::iterator it = thedron1.begin();
         it != thedron1.end(); it++) {
        int koda1 = it->first;
        Oglisca og1 = it->second;

        pair<multimap<long int, Oglisca>::iterator,
             multimap<long int, Oglisca>::iterator>
            ret;

        ret = thedron2.equal_range(koda1);
        for (multimap<long int, Oglisca>::iterator it2 = ret.first;
             it2 != ret.second; it2++) {
            Oglisca og2 = it2->second;

            /* bolj natancno primerjamo po dva tetraedra naenkrat */
            if (fabs(og1.a - og2.a) < OGA &&  // 1
                fabs(og1.b - og2.b) < OGB &&  // 1
                fabs(og1.c - og2.c) < OGC &&  // 1
                fabs(og1.V - og2.V) < OGV &&  // 2
                s->score_blosum(one_letter_code(og1.d1->atom->resn),
                                one_letter_code(og2.d1->atom->resn)) > 0 &&
                s->score_blosum(one_letter_code(og1.d2->atom->resn),
                                one_letter_code(og2.d2->atom->resn)) > 0 &&
                s->score_blosum(one_letter_code(og1.d3->atom->resn),
                                one_letter_code(og2.d3->atom->resn)) > 0 &&
                s->score_blosum(one_letter_code(og1.d4->atom->resn),
                                one_letter_code(og2.d4->atom->resn)) > 0) {
                Item* item = new Item();

                item->vert.push_back(new Vert());
                item->vert.push_back(new Vert());
                item->vert.push_back(new Vert());
                item->vert.push_back(new Vert());

                item->vert[0]->row = new Row();
                item->vert[1]->row = new Row();
                item->vert[2]->row = new Row();
                item->vert[3]->row = new Row();

                item->vert[0]->row->desc1 = og1.d1;
                item->vert[0]->row->desc2 = og2.d1;

                item->vert[1]->row->desc1 = og1.d2;
                item->vert[1]->row->desc2 = og2.d2;

                item->vert[2]->row->desc1 = og1.d3;
                item->vert[2]->row->desc2 = og2.d3;

                item->vert[3]->row->desc1 = og1.d4;
                item->vert[3]->row->desc2 = og2.d4;

                //        dbgmsg("-");
                //        int janez = 0;

                float rmsd = k->superimpose_calpha(item);

                //        cout << "tetraeder "<<endl
                //             <<og1.d1->atom->resn<<og1.d1->atom->resi<< "-"
                //             <<og2.d1->atom->resn<<og2.d1->atom->resi<<endl
                //             <<og1.d2->atom->resn<<og1.d2->atom->resi<< "-"
                //             <<og2.d2->atom->resn<<og2.d2->atom->resi<<endl
                //             <<og1.d3->atom->resn<<og1.d3->atom->resi<< "-"
                //             <<og2.d3->atom->resn<<og2.d3->atom->resi<<endl
                //             <<og1.d4->atom->resn<<og1.d4->atom->resi<< "-"
                //             <<og2.d4->atom->resn<<og2.d4->atom->resi<<endl
                //             <<og1.d1->s->mnsp<< "-" <<og2.d1->s->mnsp<<endl
                //             <<og1.d2->s->mnsp<< "-" <<og2.d2->s->mnsp<<endl
                //             <<og1.d3->s->mnsp<< "-" <<og2.d3->s->mnsp<<endl
                //             <<og1.d4->s->mnsp<< "-" <<og2.d4->s->mnsp<<endl
                //             << "rmsd_calpha = " << rmsd << endl;

                if (rmsd >= 0 && rmsd < THRMSD) {  // 1.2

                    //          cout <<  "MATRIX> " << (float)
                    //          gsl_matrix_get(item->U, 0, 0) << " ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 0, 1) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 0, 2) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 1, 0) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 1, 1) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 1, 2) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 2, 0) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 2, 1) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_matrix_get(item->U, 2, 2) << "
                    //                        ";
                    //                        cout <<  (float)
                    //                        gsl_vector_get(item->t, 0) << " ";
                    //                        cout <<  (float)
                    //                        gsl_vector_get(item->t, 1) << " ";
                    //                        dbgmsg( (float)
                    //                        gsl_vector_get(item->t, 2));

                    //          cout << "tetraeder "<<endl
                    //               <<og1.d1->atom->resn<<og1.d1->atom->resi<<
                    //               "-"
                    //               <<og2.d1->atom->resn<<og2.d1->atom->resi<<endl
                    //               <<og1.d2->atom->resn<<og1.d2->atom->resi<<
                    //               "-"
                    //               <<og2.d2->atom->resn<<og2.d2->atom->resi<<endl
                    //               <<og1.d3->atom->resn<<og1.d3->atom->resi<<
                    //               "-"
                    //               <<og2.d3->atom->resn<<og2.d3->atom->resi<<endl
                    //               <<og1.d4->atom->resn<<og1.d4->atom->resi<<
                    //               "-"
                    //               <<og2.d4->atom->resn<<og2.d4->atom->resi<<endl
                    //               <<og1.d1->s->mnsp<< "-"
                    //               <<og2.d1->s->mnsp<<endl
                    //               <<og1.d2->s->mnsp<< "-"
                    //               <<og2.d2->s->mnsp<<endl
                    //               <<og1.d3->s->mnsp<< "-"
                    //               <<og2.d3->s->mnsp<<endl
                    //               <<og1.d4->s->mnsp<< "-"
                    //               <<og2.d4->s->mnsp<<endl
                    //               << "rmsd_calpha = " << rmsd << endl;

                    i++;
                    //          cout << "tetraeder " << " d1 = " <<
                    //          og1.d1->atom->resi<< " d2 = " <<
                    //          og2.d1->atom->resi
                    //               << " d3 = " << og1.d2->atom->resi<< " d4 =
                    //               " << og2.d2->atom->resi
                    //               << " d5 = " << og1.d3->atom->resi<< " d6 =
                    //               " << og2.d3->atom->resi
                    //               << " d7 = " << og1.d4->atom->resi<< " d8 =
                    //               " << og2.d4->atom->resi
                    //               << " rmsd_calpha = " << rmsd << endl;

                    filtrat.push_back(make_pair(make_pair(og1, og2), item));
                } else {
                    delete item->vert[0]->row;
                    delete item->vert[1]->row;
                    delete item->vert[2]->row;
                    delete item->vert[3]->row;
                    delete item;
                }
            }
        }
    }

    int j = 0;
    /* izlocimo tiste tetraedre, ki dajo enako translacijo in rotacijo */
    for (vector<pair<pair<Oglisca, Oglisca>, Item*> >::iterator it =
             filtrat.begin();
         it != filtrat.end(); it++) {
        Item* mat1 = it->second;

        for (vector<pair<pair<Oglisca, Oglisca>, Item*> >::iterator it2 =
                 it + 1;
             it2 != filtrat.end(); it2++) {
            Item* mat2 = it2->second;

            if (k->compare_matrices(mat1->U, mat2->U, mat1->t, mat2->t)) {
                delete mat2->vert[0]->row;
                delete mat2->vert[1]->row;
                delete mat2->vert[2]->row;
                delete mat2->vert[3]->row;

                delete mat2;

                filtrat.erase(it2);
                it2--;
                j++;
            }
        }
    }

    delete k;

    dbgmsg("Product::tetrahedron Time = " << (double)(clock() - start) /
                                                 CLOCKS_PER_SEC);
    dbgmsg("Product::tetrahedron Number of equal tetrahedrons = " << i);
    dbgmsg("Product::tetrahedron Number of same matrices = " << j);
    dbgmsg("Product::tetrahedron Number of all possible equal tetrahedrons = "
           << thedron1.size() * thedron2.size());
}

void Product::init_product(Score* s) {
    /*
       Tukaj inicializiramo product graph.
    */

    dbgmsg("Initializing product graph ...");

    Row* r;
#ifndef NDEBUG
    clock_t start = clock();
#endif
    row = new Row();

    r = row;

    size = 0;

    Kabsch* k = new Kabsch();

    /* gremo po filtratu */
    for (vector<pair<pair<Oglisca, Oglisca>, Item*> >::iterator it =
             filtrat.begin();
         it != filtrat.end(); it++) {
        Oglisca og1 = it->first.first;
        Oglisca og2 = it->first.second;
        Item* mat = it->second;

#ifndef NDEBUG
        cout << "tetraeder[ " << size << "] d1 = " << og1.d1->atom->resi
             << " d2 = " << og2.d1->atom->resi << " d3 = " << og1.d2->atom->resi
             << " d4 = " << og2.d2->atom->resi << " d5 = " << og1.d3->atom->resi
             << " d6 = " << og2.d3->atom->resi << " d7 = " << og1.d4->atom->resi
             << " d8 = " << og2.d4->atom->resi << endl;
#endif

        // v product graf dodamo samo en (prvi) ujemajoci par deskriptorjev
        r->desc1 = og1.d1;
        r->desc2 = og2.d1;

        //    dbgmsg("DN> " << r->desc1->num << " " << r->desc2->num);

        r->num = size;

        r->debug = 0;

        /* dodamo edges */
        r->column = new Column();
        Column* c = r->column;

        int i = 1;

#ifndef NDEBUG
        Coor init_crd = k->rotate_vector(og1.d1->crd, mat->U, mat->t);
        Coor crd2 = k->rotate_vector(og1.d2->crd, mat->U, mat->t);
        Coor crd3 = k->rotate_vector(og1.d3->crd, mat->U, mat->t);
        Coor crd4 = k->rotate_vector(og1.d4->crd, mat->U, mat->t);
        dbgmsg("Razdalja med prvima deskriptorjema = " << dist(init_crd,
                                                               r->desc2->crd));
        dbgmsg(
            "Razdalja med drugim deskriptorjema = " << dist(crd2, og2.d2->crd));
        dbgmsg(
            "Razdalja med tretji deskriptorjema = " << dist(crd3, og2.d3->crd));
        dbgmsg(
            "Razdalja med cetrti deskriptorjema = " << dist(crd4, og2.d4->crd));
        cout << endl;
#endif

        /* poiscemo ujemajoce deskriptorje */
        for (Element* sosed1 = r->desc1->neighb; sosed1 != NULL;
             sosed1 = sosed1->next) {
            /* prilegamo sosede od desc1 na sosede od desc2, kakor narekuje
             * rotacijska matrika za ta filtrat */
            Descriptor* d1 = (Descriptor*)sosed1->item;

            //      dbgmsg("BB> "<< d1->bb );

            Coor neighb_crd = k->rotate_vector((d1)->crd, mat->U, mat->t);

            for (Element* sosed2 = r->desc2->neighb; sosed2 != NULL;
                 sosed2 = sosed2->next) {
                Descriptor* d2 = (Descriptor*)sosed2->item;

                //        dbgmsg("BB> "<< d2->bb );

                /* ce se ujemata v mnsp -- ce sta oba deskriptorja iz peptidne
                 * vezi, potem preverimo, ce je blosum score > 0 */
                if (d1->s->compare_single(d2->s) &&
                    (!d1->bb || !d2->bb ||
                     s->score_blosum(one_letter_code(d1->atom->resn),
                                     one_letter_code(d2->atom->resn)) > 0)) {
                    //        if (d1->s->compare_single( d2->s ) ) {

                    float dist0 = dist(neighb_crd, d2->crd);

#ifndef NDEBUG
                    if (d1->atom->resi == 844 && d2->atom->resi == 176) {
                        dbgmsg("ph1=" << d1->s->mnsp << " ph2=" << d2->s->mnsp
                                      << " at1=" << d1->atom->tag << " at2="
                                      << d2->atom->tag << " dist=" << dist0);
                    }
#endif

                    //          if (pow(dist1 + dist2 - RESOLUTION_2, 2) <
                    //          4*dist1*dist2) {
                    if (dist0 < RESOLUTION) {
                        c->item = new Row();

                        c->item->desc1 = d1;
                        c->item->desc2 = d2;

                        c->num = i++;
                        c->next = new Column();
                        c = c->next;
                    }
                }
            }
        }

        /* in narediti novi subgraph */
        r->next = new Row();
        r = r->next;
        size++;
    }

    delete k;

    dbgmsg("Product::init_product() Size of product graph: " << size);
    dbgmsg("Product::init_product() Time = " << (double)(clock() - start) /
                                                    CLOCKS_PER_SEC);
}

//    Product::Product(Descriptor *desc1, Descriptor *desc2, int psurf) {
//    //    Here we construct the product graph of the two proteins descriptors
//    sets.
//    //    First, we compare all against all vertices in the two descriptors
//    sets. If two vertices
//    //    share common neighbourhoods (>MNSP_LOW), and if both descriptors are
//    on the surface with color = lcolor,
//    //    then we create a new Row of the product graph. Each row contains the
//    two pointers to the original
//    //    descriptors and the mnsp score.
//    //    Note: psurf = 0 when total proteins search, psurf = 1 when
//    protein-protein interaction search
//      dbgmsg("Initializing product graph ...");
//      Descriptor *d1, *d2;
//      Row *r;
//      double tmp;
//      //#ifndef NDEBUG
//      int ana_score[1000], tot_score = 0;
//      int temp_score = 0;
//
//      for (int i = 0; i < 1000; i++) ana_score[i] = 0;
//
//      clock_t start = clock();
//
//      //#endif
//
//      row = new Row();
//      r = row;
//      size = 0;
//      d1 = desc1;
//      /* zacasno */
//    //  Descriptor *tmpd = desc1;
//    //  while (tmpd != NULL) {
//    //    tmpd->psurf = 0;
//    //    tmpd = tmpd->next;
//    //  }
//      /***********/
//      while (d1 != NULL) {
//        if (d1->psurf == psurf) {
//          d2 = desc2;
//          while (d2 != NULL) {
//            if (d1->s->compare_single(d2->s)) {
//              if ((tmp = d1->s->compare(d2->s)) > MNSP_LOW) {
//
//                bool dodaj = true;
//
//                /* preverimo, ali ni preblizu kateremu drugemu ze dodanemu
//                entryu product grapha */
//                for (Element *tmp2=d1->neighb; tmp2!=NULL; tmp2=tmp2->next) {
//                  float dist1 = tmp2->dAB;
//                  Row *r2;
//                  while ((r2 = ((Descriptor*) tmp2->item)->get()) != NULL) {
//                    float dist2 = dist_fast(d2->crd, r2->desc2->crd);
//                    if ((sqrt(dist1) + sqrt(dist2))/2 < 3) {
//                      dodaj = false;
//                      ((Descriptor*) tmp2->item)->reset_current();
//                      goto vun;
//                    }
//                  }
//                }
//              vun:
//                if (dodaj) {
//                  r->desc1 = d1;
//                  r->desc2 = d2;
//                  r->score = tmp;
//                  r->num = size;
//
//                  r->debug = 0;
//
//                  d1->insert(r);
//                  r->next = new Row();
//                  r = r->next;
//                  size++;
//                }
//
//              }
//    //#ifndef NDEBUG
//              ana_score[(int) (10*tmp)]++;
//              tot_score++;
//    //#endif
//              //        output_descriptor(d1, d2, tmp);
//            }
//            d2 = d2->next;
//          }
//        }
//        d1 = d1->next;
//      }
//    //#ifndef NDEBUG
//      for (int i = 0; i < 100; i++) { temp_score += ana_score[i];dbgmsg(i<<":
//      "<<ana_score[i]<<" percent:
//      "<<(100*(double)temp_score/(double)tot_score)); }
//      dbgmsg("Product::Product Total number of equ.: " << tot_score);
//      dbgmsg("Product::Product Size of product graph: " << size);
//      dbgmsg("Product::Product Time = " << (double)(clock() - start) /
//      CLOCKS_PER_SEC);
//    //#endif
//      //  exit(0);
//    }
//
//    void Product::sort_product() {
//      /*
//        Sortiramo od najvisjega do najnizjega score-a
//       */
//
//
//      Row *prev_r = NULL;
//      Row *r = row;
//      while (r->next != NULL) {
//
//        Row *prev_r1 = r;
//        Row *r1 = r->next;
//        while (r1->next != NULL) {
//
//          if (r->score < r1->score) {
//            if (prev_r != NULL) {
//              prev_r->next = r1;
//
//              Row *tmp = r1->next;
//
//              r1->next = r->next;
//
//              prev_r1->next = r;
//              r->next = tmp;
//            }
//            prev_r1 = r;
//            r1 = r->next;
//          }
//          prev_r1 = r1;
//          r1 = r1->next;
//        }
//
//        prev_r = r;
//        r = r->next;
//      }
//
//    }
//
//    void Product::init_edges() {
//      /*
//        We connect vertices with r->score score, obtained by comparison of
//        schmitt descriptors, > MNSP_HIGH.
//        The two vertices (n1,n2), (m1,m2) are connected if distance n1-m1 <
//        CUTOFF_FIRST or distance n2-m2 < CUTOFF_FIRST and
//        both distances differ maximum by RESOLUTION.
//        NOTE: Although it was not used here before this date as of JUL/27/2008
//        the CUTOFF_FIRST constant does not exist.
//      */
//      dbgmsg("Edges of the product graph are being constructed...");
//      Row *r, *r2;
//      Column *c;
//      double dist1, dist2, tmp_dist1, tmp_dist2;
//      int i;
//
//    //#ifndef NDEBUG
//      clock_t start = clock();
//    //#endif
//
//      r = row;
//      while (r->next != NULL) {
//        r->column = new Column();
//
//        dbgmsg(r->score);
//
//        if (r->score > MNSP_HIGH) {
//          i = 1;
//          c = r->column;
//          for (Element *tmp2=r->desc1->neighb; tmp2!=NULL; tmp2=tmp2->next) {
//            dist1 = tmp2->dAB;
//            while ((r2 = ((Descriptor*) tmp2->item)->get()) != NULL) {
//              dist2 = dist_fast(r->desc2->crd, r2->desc2->crd);
//              if (pow(dist1 + dist2 - RESOLUTION_2, 2) < 4*dist1*dist2) {
//                //        dbgmsg("dist1 = " << sqrt(dist1) << "dist2 = " <<
//                sqrt(dist2));
//                //        if (pow(dist1 + dist2 - resolution(dist1), 2) <
//                4*dist1*dist2) {
//                //          if (r->score > MNSP_HIGH) r2->debug = 1;
//                tmp_dist1 = sqrt(dist1);
//                tmp_dist2 = sqrt(dist2);
//                c->dist = (tmp_dist1 + tmp_dist2) / 2;
//                c->diff = fabs(tmp_dist1 - tmp_dist2);
//                c->item = r2;
//                c->num = i++;
//                c->next = new Column();
//                c = c->next;
//              }
//            }
//          }
//        }
//        r = r->next;
//      }
//      //#ifndef NDEBUG
//      dbgmsg("Product::init_edges Time = " << (double)(clock() - start) /
//      CLOCKS_PER_SEC);
//      //#endif
//      //  exit(0);
//    }

//    Product::Product(Descriptor *desc1, Descriptor *desc2, int psurf) {
//    //    Here we construct the product graph of the two proteins descriptors
//    sets.
//    //    First, we compare all against all vertices in the two descriptors
//    sets. If two vertices
//    //    share common neighbourhoods (>MNSP_LOW), and if both descriptors are
//    on the surface with color = lcolor,
//    //    then we create a new Row of the product graph. Each row contains the
//    two pointers to the original
//    //    descriptors and the mnsp score.
//    //    Note: psurf = 0 when total proteins search, psurf = 1 when
//    protein-protein interaction search
//      dbgmsg("Initializing product graph ...");
//      Descriptor *d1, *d2;
//      Row *r;
//      double tmp;
//      //#ifndef NDEBUG
//      int ana_score[1000], tot_score = 0;
//      int temp_score = 0;
//
//      for (int i = 0; i < 1000; i++) ana_score[i] = 0;
//
//      clock_t start = clock();
//
//      //#endif
//
//      row = new Row();
//      r = row;
//      size = 0;
//      d1 = desc1;
//      /* zacasno */
//    //  Descriptor *tmpd = desc1;
//    //  while (tmpd != NULL) {
//    //    tmpd->psurf = 0;
//    //    tmpd = tmpd->next;
//    //  }
//      /***********/
//      while (d1 != NULL) {
//        if (d1->psurf == psurf) {
//          d2 = desc2;
//          while (d2 != NULL) {
//            if (d1->s->compare_single(d2->s)) {
//              if ((tmp = d1->s->compare(d2->s)) > MNSP_LOW) {
//                r->desc1 = d1;
//                r->desc2 = d2;
//                r->score = tmp;
//                r->num = size;
//
//                r->debug = 0;
//
//                d1->insert(r);
//                r->next = new Row();
//                r = r->next;
//                size++;
//              }
//    //#ifndef NDEBUG
//              ana_score[(int) (10*tmp)]++;
//              tot_score++;
//    //#endif
//              //        output_descriptor(d1, d2, tmp);
//            }
//            d2 = d2->next;
//          }
//        }
//        d1 = d1->next;
//      }
//    //#ifndef NDEBUG
//      for (int i = 0; i < 100; i++) { temp_score += ana_score[i];dbgmsg(i<<":
//      "<<ana_score[i]<<" percent:
//      "<<(100*(double)temp_score/(double)tot_score)); }
//      dbgmsg("Product::Product Total number of equ.: " << tot_score);
//      dbgmsg("Product::Product Size of product graph: " << size);
//      dbgmsg("Product::Product Time = " << (double)(clock() - start) /
//      CLOCKS_PER_SEC);
//    //#endif
//      //  exit(0);
//    }
//
//    void Product::init_edges() {
//      /*
//        We connect vertices with r->score score, obtained by comparison of
//        schmitt descriptors, > MNSP_HIGH.
//        The two vertices (n1,n2), (m1,m2) are connected if distance n1-m1 <
//        CUTOFF_FIRST or distance n2-m2 < CUTOFF_FIRST and
//        both distances differ maximum by RESOLUTION.
//        NOTE: Although it was not used here before this date as of JUL/27/2008
//        the CUTOFF_FIRST constant does not exist.
//      */
//      dbgmsg("Edges of the product graph are being constructed...");
//      Row *r, *r2;
//      Column *c;
//      double dist1, dist2, tmp_dist1, tmp_dist2;
//      int i;
//
//    //#ifndef NDEBUG
//      clock_t start = clock();
//    //#endif
//
//      r = row;
//      while (r->next != NULL) {
//        r->column = new Column();
//        if (r->score > MNSP_HIGH) {
//          i = 1;
//          c = r->column;
//          for (Element *tmp2=r->desc1->neighb; tmp2!=NULL; tmp2=tmp2->next) {
//            dist1 = tmp2->dAB;
//            while ((r2 = ((Descriptor*) tmp2->item)->get()) != NULL) {
//              dist2 = dist_fast(r->desc2->crd, r2->desc2->crd);
//              if (pow(dist1 + dist2 - RESOLUTION_2, 2) < 4*dist1*dist2) {
//                //        dbgmsg("dist1 = " << sqrt(dist1) << "dist2 = " <<
//                sqrt(dist2));
//                //        if (pow(dist1 + dist2 - resolution(dist1), 2) <
//                4*dist1*dist2) {
//                //          if (r->score > MNSP_HIGH) r2->debug = 1;
//                tmp_dist1 = sqrt(dist1);
//                tmp_dist2 = sqrt(dist2);
//                c->dist = (tmp_dist1 + tmp_dist2) / 2;
//                c->diff = fabs(tmp_dist1 - tmp_dist2);
//                c->item = r2;
//                c->num = i++;
//                c->next = new Column();
//                c = c->next;
//              }
//            }
//          }
//        }
//        r = r->next;
//      }
//      //#ifndef NDEBUG
//      dbgmsg("Product::init_edges Time = " << (double)(clock() - start) /
//      CLOCKS_PER_SEC);
//      //#endif
//      //  exit(0);
//    }

//    void Product::combined_score() {
//      Row *r = row;
//      Column *c;
//      double score;
//      int i;
//      while (r->next != NULL) {
//        c = r->column;
//        score = 0.0;
//        i = 0;
//        while (c->next != NULL) {
//          score += c->item->score;
//          i++;
//          c = c->next;
//        }
//        r->combined_score = score/i;
//        r = r->next;
//      }
//    //#ifdef DEBUG
//    //  deb_out_combined_score(row);
//    //#endif
//    }

void Product::count_subgraphs() {
    Row* r = row;
    num_subgraph = 0;
    while (r->next != NULL) {
        //    if (r->score > MNSP_HIGH) {
        r->sgraph = num_subgraph;
        num_subgraph++;
        //    }
        r = r->next;
    }
    dbgmsg("Number of subgraphs = " << num_subgraph);
}

// void Product::remove_contacts() {
//  Row *r = row;
//  Column *c;
//  while (r->next != NULL) {
//    if (r->score > MNSP_HIGH) {
//      c = r->column;
//      while (c->next != NULL) {
//        if (c->dist < CUTOFF_SECOND)
//          if (c->item->score > MNSP_HIGH) {
//            c->item->score = -(c->item->score);
//          }
//        c = c->next;
//      }
//    }
//    r = r->next;
//  }
//  r = row;
//  while (r->next) {
//    if (r->score > MNSP_HIGH) {
//      c = r->column;
//      r->debug = 1;
//      while (c->next) {
//        c->item->debug = 1;
//        c = c->next;
//      }
//    }
//    r = r->next;
//  }
//  r = row;
//  int a = 0, b = 0;
//  while (r->next) {
//    if (r->debug)
//      a++;
//    else
//      b++;
//    r = r->next;
//  }
//  dbgmsg("number of rows with debug = 1 : " << a << " and rows with debug = 0
//  : " << b);
////  exit(0);
//}

// void Product::remove_helix() {
//  dbgmsg("Remove alpha helices...");
//  Row *r = row;
//  while (r->next != NULL) {
//    if (r->score > MNSP_HIGH) {
//      if (r->desc1->atom->helix && r->desc2->atom->helix)
//          //      if (((r->desc1->atom->helix && backbone) ||
//          (r->desc2->atom->helix && backbone))
//        r->score = - (r->score);
//    }
//    r = r->next;
//  }
//}

//    void Product::tetraedri() {
//      /*
//        Poiscemo podobne tetraedre v obeh proteinih, pri cemer so oglisca
//        deskriptorji z enakimi oznakami.
//      */
//      dbgmsg("Edges of the product graph are being constructed...");
//
//    //#ifndef NDEBUG
//      clock_t start = clock();
//    //#endif
//
//      Row *r1 = row;
//      int i = 0;
//
//      dbgmsg("RESOLUTION_2 = " << RESOLUTION_2);
//
//      while (r1->next != NULL) {
//        //    r1->column = new Column();
//
//        dbgmsg("i =  " << i++);
//
//        for (Element *d2=r1->desc1->neighb; d2!=NULL; d2=d2->next) {
//
//          double dist1 = d2->dAB;
//          Row *r2;
//
//          while ((r2 = ((Descriptor*) d2->item)->get()) != NULL) {
//    //        dbgmsg("prva zanka ");
//
//              double dist2 = dist_fast(r1->desc2->crd, r2->desc2->crd);
//
//              if (pow(dist1 + dist2 - RESOLUTION_2, 2) < 4*dist1*dist2) {
//
//                for (Element *d3=d2->next; d3!=NULL; d3=d3->next) {
//
//                  double dist3 = d3->dAB;
//
//                  Row *r3;
//
//                  while ((r3 = ((Descriptor*) d3->item)->get()) != NULL) {
//    //                dbgmsg("druga zanka ");
//
//                    double dist4 = dist_fast(r1->desc2->crd, r3->desc2->crd);
//
//                    if (pow(dist3 + dist4 - RESOLUTION_2, 2) < 4*dist3*dist4)
//                    {
//
//                      double dist5 = dist_fast(r2->desc1->crd,
//                      r3->desc1->crd);
//                      double dist6 = dist_fast(r2->desc2->crd,
//                      r3->desc2->crd);
//
//                      if (pow(dist5 + dist6 - RESOLUTION_2, 2) <
//                      4*dist5*dist6) {
//
//
//                        for (Element *d4=d3->next; d4!=NULL; d4=d4->next) {
//
//                          double dist7 = d4->dAB;
//
//                          Row *r4;
//
//                          while ((r4 = ((Descriptor*) d4->item)->get()) !=
//                          NULL) {
//                            //                      dbgmsg("tretja zanka ");
//
//                          double dist8 = dist_fast(r1->desc2->crd,
//                          r4->desc2->crd);
//
//                          if (pow(dist7 + dist8 - RESOLUTION_2, 2) <
//                          4*dist7*dist8) {
//
//                            double dist9  = dist_fast(r2->desc1->crd,
//                            r4->desc1->crd);
//                            double dist10 = dist_fast(r2->desc2->crd,
//                            r4->desc2->crd);
//
//                            if (pow(dist9 + dist10 - RESOLUTION_2, 2) <
//                            4*dist9*dist10) {
//
//                              double dist11 = dist_fast(r3->desc1->crd,
//                              r4->desc1->crd);
//                              double dist12 = dist_fast(r3->desc2->crd,
//                              r4->desc2->crd);
//
//                              if (pow(dist11 + dist12 - RESOLUTION_2, 2) <
//                              4*dist11*dist12) {
//                                cout << "nasel tetraeder "
//                                     << " d1 = " << sqrt(dist1) << " d2 = " <<
//                                     sqrt(dist2)
//                                     << " d3 = " << sqrt(dist3) << " d4 = " <<
//                                     sqrt(dist4)
//                                     << " d5 = " << sqrt(dist5) << " d6 = " <<
//                                     sqrt(dist6)
//                                     << " d7 = " << sqrt(dist7) << " d8 = " <<
//                                     sqrt(dist8)
//                                     << " d9 = " << sqrt(dist9) << " d10 = "<<
//                                     sqrt(dist10)
//                                     << " d11 = " << sqrt(dist11) << " d12 = "
//                                     << sqrt(dist12)
//                                     << endl;
//                                cout << "nasel tetraeder "
//                                     << " d1 = " << r1->desc1->atom->resi<< "
//                                     d2 = " << r1->desc2->atom->resi
//                                     << " d3 = " << r2->desc1->atom->resi<< "
//                                     d4 = " << r2->desc2->atom->resi
//                                     << " d5 = " << r3->desc1->atom->resi<< "
//                                     d6 = " << r3->desc2->atom->resi
//                                     << " d7 = " << r4->desc1->atom->resi<< "
//                                     d8 = " << r4->desc2->atom->resi
//                                     << endl;
//                              }
//                            }
//                          }
//                          }
//                        }
//                      }
//                    }
//                  }
//                }
//              }
//          }
//        }
//        r1 = r1->next;
//      }
//
//      dbgmsg("Product::tetraedri Time = " << (double)(clock() - start) /
//      CLOCKS_PER_SEC);
//
//    }
