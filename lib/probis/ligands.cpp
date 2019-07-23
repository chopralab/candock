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

#include "ligands.h"
#include "bsite.h"
#include "clusterdata.h"
#include "clusterlib.h"
#include "geo.h"
#include "grid.h"
#include "kabsch.h"
#include "molecule.h"

Ligands::~Ligands() {
    /* sprostimo pomnilnik za ligands */
    for (multimap<Residue*, Ligand*>::iterator it = ligands.begin();
         it != ligands.end(); it++) {
        delete it->first;
        delete it->second;
    }

    /* sprostimo pomnilnik za lig_query */
    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        delete *it;
    }

    // sprostimo data
    for (int i = 0; i < nPP; i++) delete[] dataPP[i];
    for (int i = 0; i < nSL; i++) delete[] dataSL[i];
    for (int i = 0; i < nNU; i++) delete[] dataNU[i];
    for (int i = 0; i < nIO; i++) delete[] dataIO[i];
    delete[] dataPP;
    delete[] dataSL;
    delete[] dataNU;
    delete[] dataIO;

    // sprostimo plig
    delete[] pligPP;
    delete[] pligSL;
    delete[] pligNU;
    delete[] pligIO;

    // sprostimo zazipani data in plig (same vsebine zdata ni treba sprostit,
    // ker jo ze pri data)
    delete[] zdataPP;
    delete[] zdataSL;
    delete[] zdataNU;
    delete[] zdataIO;

    delete[] zpligPP;
    delete[] zpligSL;
    delete[] zpligNU;
    delete[] zpligIO;
}

void Ligands::read(Molecule* m) {
    /*
       Iz datotek v ustreznem direktoriju LIGDIR preberemo ligande, ki smo jih
       nasli v vsakem klastru.
       Input datoteke smo dobili tako, da smo ligande iz vseh PDB-jev v klastru
       prilegali na predstojnika
       (PDB z najvecjo locljivostjo) in za vsak ligand smo nasli njemu bliznje
       (<3A) aminokisline v predstojniku.

       Output je multimap tabela ligands oblike:

       PREDSTOJNIK         LIGAND
       -------------------------------------
       pdb_id,resi         pdb_id,resn,resi,chain
       -------------------------------------
       1i52A,227     -->   1injA,HOH,1033,A
       .                   .
       .                   .
       .                   .

    */

    /* gremo po vseh (najvec prvih NMOL) prileganih proteinih
     */
    int i = 0;
    for (vector<Molecule*>::iterator it = m->lista_molekul.begin();
         i < NMOL && it != m->lista_molekul.end(); it++) {
        i++;

        //    string align_pdb = to_string((*it)->pdb_id);
        string align_pdb =
            to_string((*it)->pdb_id) + to_string((*it)->chain_id);

#ifdef VERB
        cout << "Ligands::read align_pdb = " << align_pdb << endl;
#endif

        /* preberemo ustrezno datoteko tipa .lig (npr. 1i52A.lig), formata:

           PREDSTOJNIK                                            LIGAND IZ
           PRILEGANEGA
           ----------------------------------------------------------------------------------------------------
                #MODEL     PDB        RESN     #RESI     CHAIN    #MODEL     PDB
           RESN     #RESI      CHAIN
           ----------------------------------------------------------------------------------------------------
           01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
           LIG>      1     1i52A       LEU       227         A         5
           1injA       HOH      1033         A
           LIG>      1     1i52A       LEU       227         A         5
           1injA       @@@         1         A
           ----------------------------------------------------------------------------------------------------

           kadar je ligand protein, ima oznake resn=@@@ in resi=1 :)
           kadar je ligand nucleic, ima oznake resn=+++ in resi=1 :)

        */

        string name = add_end_slash(LIGDIR) + "ligands/" + align_pdb + ".lig";

        ifstream f(name.c_str());
        string line;

        if (f.is_open()) {
            while (!f.eof()) {
                line.clear();

                getline(f, line);

                if (line.substr(0, 4).compare("LIG>") == 0) {
                    //          string pdb_id1 = line.substr(16, 5);
                    int resi1 = atoi(line.substr(36, 5).c_str());
                    char chain_id1 = line.at(50);

                    /* Na primer, [model=3,1xxyA,ATP,301,B] in
                     * [model=4,1xxyA,ATP,301,B] sta obravnavana kot dva
                     * razlicna liganda */
                    int model2 = atoi(line.substr(55, 6).c_str());
                    //         string pdb_id2 = line.substr(66, 5);
                    string pdb_id2 = line.substr(66, 4);
                    char acid =
                        line.at(70);  // NOVO  // chain_id alignane verige
                    string resn2 = line.substr(78, 3);
                    int resi2 = atoi(line.substr(86, 5).c_str());
                    char chain_id2 = line.at(100);

                    // zapisemo tudi tip liganda
                    bstype t;
                    if (resn2.compare("@@@") == 0) {
                        t = _pp;
                    } else if (resn2.compare("+++") == 0) {
                        t = _nu;
                    } else if (ions.find(resn2) != string::npos) {
                        t = _io;
                    } else {
                        t = _sl;
                    }

                    /* vode in site zapisov (resi2 = 0) ne upostevamo
                     * potencialna napaka : ligand ima tudi lahko resi = 0 -
                     * ceprav precej neverjetno :-) */
                    if (resn2.compare("HOH") != 0 &&
                        !(resi2 == 0 ||
                          resi2 == -999)) {  // po novem je SITE zapis = -999 (
                                             // zaradi kompatibilnosti imamo se
                                             // vedno resi2 == 0 )
                        //	  if (resn2.compare("HOH") != 0 &&
                        //!(resn2.compare(0,2, "AC") == 0 && (resi2 == 0 ||
                        //resi2 == -999))) {
                        /* in naredimo tabelo ligands */
                        //            ligands.insert( make_pair( new Residue(
                        //            chain_id1, resi1, *it ), new Ligand(
                        //            model2, pdb_id2, chain_id2, resn2, resi2 ,
                        //            t) ) );
                        ligands.insert(make_pair(
                            new Residue(chain_id1, resi1, *it),
                            new Ligand(model2, pdb_id2, acid, chain_id2, resn2,
                                       resi2, t)));  // NOVO
                    }
                }
            }

        } else {
            cout << "Warning (MOLECULE) : Missing ligand file " << name << endl;
        }

        f.close();
    }
#ifdef VERB
    for (multimap<Residue*, Ligand*>::iterator it = ligands.begin();
         it != ligands.end(); it++) {
        /* s2 predstavlja ligand in je edinstven za vsak ligand */
        Residue* r1 = it->first;
        Ligand* r2 = it->second;

        cout << "Ligands::read() ligands " << r1->m->pdb_id << " "
             << r1->m->chain_id << " " << r1->resi << " " << r2->pdb_id << " "
             << r2->resn << " " << r2->resi << endl;
    }
#endif
}

void Ligands::generate(Molecule* m) {
    /*
      Naredimo multimap tabelo (liq_query), ki vsak ligand povezuje s query
      residue-ji (lahko vecimi).
      Ligande 'poberemo' iz prileganih podobnih proteinov.

      razlicni ligandi so lahko 'vezani' na eno samo ak.

      align_pdb     align_resi                  resn
      1 1i52A   LEU   227     A -->      5 1injA   HOH  1033     A
      1 1i52A   LEU   227     A -->      7 1vgtB   HOH   271     B

    */

    //#ifdef VERB
    clock_t start = clock();
    //#endif

    /* gremo po vseh (prvih NMOL najvec) prileganih proteinih */
    int i = 0;
    for (vector<Molecule*>::iterator it = m->lista_molekul.begin();
         i < NMOL && it != m->lista_molekul.end(); it++) {
        i++;

#ifdef VERB
        cout << "Ligands::generate() align_pdb = " << (*it)->pdb_id << endl;
#endif

        // gremo po zaporedju query proteina
        for (set<Residue*>::iterator it2 = m->sequence.begin();
             it2 != m->sequence.end(); it2++) {
            /* kandidatov za mesto v query zaporedju je vec, izberemo pa tisti
               klaster z najboljsim (najvisjim) cluster_score
               in flx (=flexible) == false, ce je mozno (torej, tisto
               ne-flexible z najvisjim cluster_score) */
            pair<multimap<Residue*, Residue*>::iterator,
                 multimap<Residue*, Residue*>::iterator>
                ret = (*it)->align.equal_range(*it2);

            if (ret.first != ret.second) {
                // aligned residue, po moznosti non-fleksible in z najvisjim
                // cluster_score
                Residue* rmax =
                    max_element(ret.first, ret.second, by_noflx_cluster_score)
                        ->second;
                ChResi qres = ChResi((*it2)->chain_id, (*it2)->resi);
                //        int align_resi = rmax->resi;
                ClusterData* c = rmax->cd;

#ifdef VERB
                cout << "Ligands::generate() "
                     << "(" << qres.first << "," << qres.second << ") --> ("
                     << rmax->resi << "," << c->cluster_score << ")" << endl;
#endif

                /* gremo po vseh razlicnih ligandih za eno prilegano
                 * aminokislino (tisto iz najboljsega klastra) */
                pair<multimap<Residue*, Ligand*>::iterator,
                     multimap<Residue*, Ligand*>::iterator>
                    ret2 = ligands.equal_range(rmax);

                for (multimap<Residue*, Ligand*>::iterator it3 = ret2.first;
                     it3 != ret2.second; it3++) {
                    // ali ligand ze obstaja? (vsak razlicen cluster_id je svoj
                    // ligand)
                    it3->second->cd = c;
                    // en ligand, npr. ATP 301 lahko pripada vecim klastrom 0,
                    // 1, itd.. (saj probis najde vec razlicnih lokal
                    // alignmentov)
                    // posortirani po type in znotraj enega tipa po size_neighb
                    set<Ligand*, ligcomp>::iterator tmp_it =
                        lig_query.find(it3->second);

                    if (tmp_it == lig_query.end()) {
                        // ce liganda ni, naredimo novega
                        Ligand* li = new Ligand();
                        *li = *it3->second;  // ker je lahko en ligand v vec
                                             // razlicnih cluster_id-jih moramo
                                             // narediti kopijo liganda !
                        // novemu ligandu dodamo nov cluster_id in query residue
                        li->cd = c;
                        li->rlist.insert(qres);
                        lig_query.insert(li);
                    } else {
                        // obstojecemu ligandu dodamo nov query residue
                        (*tmp_it)->rlist.insert(qres);
                    }
                }
            }
        }
    }

    //#ifdef VERB
    cout << "Ligands::generate() Time = "
         << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
//#endif
#ifdef VERB

    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        // izpisemo rotacijsko matriko in translacijski vektor
        cout << "Ligand::generate() MATRIX "
             << gsl_matrix_get((*it)->cd->U, 0, 0) << " "
             << gsl_matrix_get((*it)->cd->U, 0, 1) << " "
             << gsl_matrix_get((*it)->cd->U, 0, 2) << " " << endl;
        cout << "Ligand::generate() MATRIX "
             << gsl_matrix_get((*it)->cd->U, 1, 0) << " "
             << gsl_matrix_get((*it)->cd->U, 1, 1) << " "
             << gsl_matrix_get((*it)->cd->U, 1, 2) << " " << endl;
        cout << "Ligand::generate() MATRIX "
             << gsl_matrix_get((*it)->cd->U, 2, 0) << " "
             << gsl_matrix_get((*it)->cd->U, 2, 1) << " "
             << gsl_matrix_get((*it)->cd->U, 2, 2) << " " << endl;

        cout << "Ligand::generate() VECTOR " << gsl_vector_get((*it)->cd->t, 0)
             << " " << gsl_vector_get((*it)->cd->t, 1) << " "
             << gsl_vector_get((*it)->cd->t, 2) << " " << endl;

        // gremo po listi residuejev, ki pripadajo enemu ligandu
        for (set<ChResi>::iterator it2 = (*it)->rlist.begin();
             it2 != (*it)->rlist.end(); it2++) {
            ChResi res = *it2;
            cout << "Ligands::generate() lig_query " << (*it)->resi << " "
                 << (*it)->pdb_id << " " << (*it)->acid << " "
                 << (*it)->chain_id << " " << (*it)->model << " " << (*it)->resn
                 << " " << (*it)->cd->cluster_id << " --> " << res.first << " "
                 << res.second << endl;
        }
    }
#endif
}

void Ligands::init_grid_ligands(Grid* grid, float distance) {
    /*
       Naredimo grid za ligande (za vsak tip liganda posebej) in poiscemo sosede
       od vsakega centra liganda.
    */

    // gremo po vseh ligandih
    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        /* v gridu vsi ligandi, ki so v mejah (0,0,0)
         * (NUM_CELLS,NUM_CELLS,NUM_CELLS) - make_grid sam preveri!!! */
        grid->make_grid(*it);
    }

    /* naredimo listo neighborjev */
    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        grid->volume_slice_lig(*it, distance, distance);
    }
}

void Ligands::center(Molecule* m) {
    /*
       Izracunamo centre (oz medoide) binding site residuejev za ligande. Na
       primer, imamo ligand ATP,
       bliznji residueji na query-ju so 156,157,158, izracunamo medoido teh treh
       residuejev.
    */

    /* preberemo plain_pdb query proteina in povezemo stevilko aminokisline s
     * koordinato njenega Calpha atoma */
    map<ChResi, Coor> resi_to_coor;

    /* gremo po CALPHA atomih query chain-a, upostevamo samo atome proteina */
    for (int i = 0; i < (int)m->plain_pdb.size(); i++) {
        if (m->plain_pdb[i].compare(0, 4, "ATOM") == 0 &&
            m->chain_id.find(m->plain_pdb[i].at(21)) != string::npos &&
            m->plain_pdb[i].compare(13, 2, "CA") == 0 &&
            amino.find(m->plain_pdb[i].substr(17, 3)) != string::npos) {
            ChResi chr = ChResi(m->plain_pdb[i].at(21),
                                atoi(m->plain_pdb[i].substr(22, 4).c_str()));

            Coor c;
            c.x = atof(m->plain_pdb[i].substr(30, 8).c_str());
            c.y = atof(m->plain_pdb[i].substr(38, 8).c_str());
            c.z = atof(m->plain_pdb[i].substr(46, 8).c_str());

            resi_to_coor[chr] = c;
        }
    }

    // gremo po ligandih, le ti so sortirani po ligcomp
    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        Coor c;
        set_zero(c);

        // gremo po listi residuejev, ki pripadajo temu ligandu
        for (set<ChResi>::iterator it2 = (*it)->rlist.begin();
             it2 != (*it)->rlist.end(); it2++) {
            ChResi qres = *it2;
            c = c + resi_to_coor[qres];

#ifdef VERB
            cout << "Ligands::center() lig_query " << (*it)->pdb_id
                 << (*it)->acid << " " << (*it)->resn << " " << (*it)->resi
                 << " " << (*it)->chain_id << " " << (*it)->cd->cluster_id
                 << " --> " << qres.first << " " << qres.second << endl;
#endif
        }
        // povprecimo koordinate
        c = c * (1 / (float)(*it)->rlist.size());
        // povezemo ligand in center njegovih binding site residuejev : ATP -->
        // (x,y,z)
        (*it)->crd = c;

#ifdef VERB
        cout << "Ligands::center() coor " << (*it)->pdb_id << " " << (*it)->acid
             << " " << (*it)->resn << " " << (*it)->resi << " "
             << (*it)->chain_id << " " << (*it)->cd->cluster_id << " --> "
             << "(" << c.x << "," << c.y << "," << c.z << ")" << endl;
#endif
    }
}

void Ligands::initialize_data() {
    /*
      Pripravimo data array za vsak tip liganda, ki hrani koordinate vseh
      ligandov, nato array plig, ki hrani pointerje na
      Ligande in povezuje koordinate v data arrayu z setom ligandov.
    */

    // posortiramo ligande po tipu in znotraj posameznega tipa po stevilu
    // sosedov
    vector<Ligand*> lig_query2;
    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        lig_query2.push_back(*it);
    }
    sort(lig_query2.begin(), lig_query2.end(), by_type_size_neighb);

    // koliko je posameznih ligandov po tipih
    for (vector<Ligand*>::iterator it = lig_query2.begin();
         it != lig_query2.end(); it++) {
        switch ((*it)->type) {
            case _pp:
                nPP++;
                break;
            case _sl:
                nSL++;
                break;
            case _nu:
                nNU++;
                break;
            case _io:
                nIO++;
                break;
        }
    }

    //#ifdef VERB
    cout << "Ligands::initialize_data() nPP = " << nPP << " nSL = " << nSL
         << " nNU = " << nNU << " nIO = " << nIO << endl;
    //#endif

    /* v data so vse koordinate vseh najdenih ligandov za vsak tip liganda
     * posebej (protein,dna/rna,small-ligand,ion)
     */
    dataPP = new double*[nPP];
    dataSL = new double*[nSL];
    dataNU = new double*[nNU];
    dataIO = new double*[nIO];

    /* v plig so pointerji nazaj na class Ligand v lig_query
     */
    pligPP = new Ligand*[nPP];
    pligSL = new Ligand*[nSL];
    pligNU = new Ligand*[nNU];
    pligIO = new Ligand*[nIO];

    /* napolnimo data array-e za treecluster (hierarchical clustering algorithm)
     */
    int iPP, iSL, iNU, iIO;
    iPP = iSL = iNU = iIO = 0;

    for (vector<Ligand*>::iterator it = lig_query2.begin();
         it != lig_query2.end(); it++) {
        Coor bsc = (*it)->crd;
        // zapisemo koordinate za vsak tip liganda
        switch ((*it)->type) {
            case _pp:
                dataPP[iPP] = new double[ncolumns];
                dataPP[iPP][0] = bsc.x;
                dataPP[iPP][1] = bsc.y;
                dataPP[iPP][2] = bsc.z;
                pligPP[iPP] = *it;
                iPP++;
                break;
            case _sl:
                dataSL[iSL] = new double[ncolumns];
                dataSL[iSL][0] = bsc.x;
                dataSL[iSL][1] = bsc.y;
                dataSL[iSL][2] = bsc.z;
                pligSL[iSL] = *it;
                iSL++;
                break;
            case _nu:
                dataNU[iNU] = new double[ncolumns];
                dataNU[iNU][0] = bsc.x;
                dataNU[iNU][1] = bsc.y;
                dataNU[iNU][2] = bsc.z;
                pligNU[iNU] = *it;
                iNU++;
                break;
            case _io:
                dataIO[iIO] = new double[ncolumns];
                dataIO[iIO][0] = bsc.x;
                dataIO[iIO][1] = bsc.y;
                dataIO[iIO][2] = bsc.z;
                pligIO[iIO] = *it;
                iIO++;
                break;
        }
    }
}

void Ligands::zip_points() {
    /*
       Zdruzimo bliznje tocke (sosede), ki so manj kot X Angstromov narazen
       vrnemo pointer na novo alociran data array,
       ki je manjsi (zazipan) vhodni (data).
    */

    zdataPP = new double*[nPP];
    zdataSL = new double*[nSL];
    zdataNU = new double*[nNU];
    zdataIO = new double*[nIO];

    zpligPP = new Ligand*[nPP];
    zpligSL = new Ligand*[nSL];
    zpligNU = new Ligand*[nNU];
    zpligIO = new Ligand*[nIO];

    // zdruzimo (greedy) sosednje tocke - vsaka tocka je lahko upostevana le
    // enkrat
    for (int i = 0; i < nPP; i++) {
        if (!pligPP[i]->visited) {
            pligPP[i]->visited = true;
            zpligPP[znPP] = pligPP[i];
            zdataPP[znPP] = dataPP[i];
#ifdef VERB
            cout << "Ligands::zip_points  :" << zpligPP[znPP]->resn << endl;
#endif
            znPP++;
#ifdef VERB
            if (!pligPP[i]->first_neighb()) {
                cout << "Ligands::zip_points no first neighbor" << endl;
                exit(0);
            }
#endif
            for (Ligand* li = pligPP[i]->first_neighb(); li != NULL;
                 li = pligPP[i]->next_neighb()) {
                if (li->visited) {
                    pligPP[i]->delete_neighb();
                } else
                    li->visited = true;
            }
        }
    }
    for (int i = 0; i < nSL; i++) {
        if (!pligSL[i]->visited) {
            pligSL[i]->visited = true;
            zpligSL[znSL] = pligSL[i];
            zdataSL[znSL] = dataSL[i];
#ifdef VERB
            cout << "Ligands::zip_points  :" << zpligSL[znSL]->resn << endl;
#endif
            znSL++;
#ifdef VERB
            if (!pligSL[i]->first_neighb()) {
                cout << "Ligands::zip_points no first neighbor" << endl;
                exit(0);
            }
#endif
            for (Ligand* li = pligSL[i]->first_neighb(); li != NULL;
                 li = pligSL[i]->next_neighb()) {
                if (li->visited) {
                    pligSL[i]->delete_neighb();
                } else
                    li->visited = true;
            }
        }
    }
    for (int i = 0; i < nNU; i++) {
        if (!pligNU[i]->visited) {
            pligNU[i]->visited = true;
            zpligNU[znNU] = pligNU[i];
            zdataNU[znNU] = dataNU[i];
#ifdef VERB
            cout << "Ligands::zip_points  :" << zpligNU[znNU]->resn << endl;
#endif
            znNU++;
#ifdef VERB
            if (!pligNU[i]->first_neighb()) {
                cout << "Ligands::zip_points no first neighbor" << endl;
                exit(0);
            }
#endif
            for (Ligand* li = pligNU[i]->first_neighb(); li != NULL;
                 li = pligNU[i]->next_neighb()) {
                if (li->visited) {
                    pligNU[i]->delete_neighb();
                } else
                    li->visited = true;
            }
        }
    }
    for (int i = 0; i < nIO; i++) {
        if (!pligIO[i]->visited) {
            pligIO[i]->visited = true;
            zpligIO[znIO] = pligIO[i];
            zdataIO[znIO] = dataIO[i];
#ifdef VERB
            cout << "Ligands::zip_points  :" << zpligIO[znIO]->resn << endl;
#endif
            znIO++;
#ifdef VERB
            if (!pligIO[i]->first_neighb()) {
                cout << "Ligands::zip_points no first neighbor" << endl;
                exit(0);
            }
#endif
            for (Ligand* li = pligIO[i]->first_neighb(); li != NULL;
                 li = pligIO[i]->next_neighb()) {
                if (li->visited) {
                    pligIO[i]->delete_neighb();
                } else
                    li->visited = true;
            }
        }
    }

#ifdef VERB
    cout << "Ligands::zip_points  : znPP = " << znPP << " znSL = " << znSL
         << " znNU = " << znNU << " znIO = " << znIO << endl;
#endif
}

void Ligands::cluster() {
    /*
       Klastriramo centre (x,y,z) ligandov, tako da bliznji v prostoru pridejo v
       isti klaster.
       Pri tem locimo klastre za proteinske ligande, za ligande, ki so
       nukleinske kisline, za
       male ligande (nekaj atomov) in enoatomne ligande (ione).
     */

    /* pripravimo za klastriranje ligandov, tako da bodo tisti, blizu v
     * prostoru, v istem klastru */
    int nrows = lig_query.size();
    double weight[3] = {1.0, 1.0, 1.0};  // [ncolumns] za transpose == 0

    /* mask je en za vse ligande (je namerno dolg, kolikor je vseh ligandov
     * skupaj)
     */
    int** mask = new int*[nrows];

    // maske inicializiramo
    for (int i = 0; i < nrows; i++) {
        mask[i] = new int[ncolumns];
        for (int j = 0; j < ncolumns; j++) mask[i][j] = 1;
    }

    /* pozenemo tree - klastering za vsak tip liganda posebej
     */
    //#ifdef VERB
    clock_t start = clock();
    //#endif
    Node* NPP =
        treecluster(znPP, ncolumns, zdataPP, mask, weight, 0, 'e', 'c', NULL);
    Node* NSL =
        treecluster(znSL, ncolumns, zdataSL, mask, weight, 0, 'e', 'c', NULL);
    Node* NNU =
        treecluster(znNU, ncolumns, zdataNU, mask, weight, 0, 'e', 'c', NULL);
    Node* NIO =
        treecluster(znIO, ncolumns, zdataIO, mask, weight, 0, 'e', 'c', NULL);
    //#ifdef VERB
    cout << "Ligands::cluster() Time = "
         << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
    cout << "Ligands::cluster() nrows = " << nrows << endl;
    //#endif

    /* skusamo dobiti optimalno stevilo klastrov, seveda za vsak tip liganda
     * posebej
     */
    int* clusidPP = new int[znPP];
    int* clusidSL = new int[znSL];
    int* clusidNU = new int[znNU];
    int* clusidIO = new int[znIO];
    int* maskidPP = new int[znPP];
    int* maskidSL = new int[znSL];
    int* maskidNU = new int[znNU];
    int* maskidIO = new int[znIO];

    /* inicializiramo (na 0) vse elemente v maskidXYZ, arrayu, ki odgovarja
       clusidXYZ, negativna vrednost pomeni, da
       klastra, ki je v clusid na tem mestu, ne delimo vec naprej.
    */
    for (int i = 0; i < znPP; i++) {
        clusidPP[i] = 0;
        maskidPP[i] = 0;
    }
    for (int i = 0; i < znSL; i++) {
        clusidSL[i] = 0;
        maskidSL[i] = 0;
    }
    for (int i = 0; i < znNU; i++) {
        clusidNU[i] = 0;
        maskidNU[i] = 0;
    }
    for (int i = 0; i < znIO; i++) {
        clusidIO[i] = 0;
        maskidIO[i] = 0;
    }

    /* zacnemo pri enem klastru (najvecjem) in postopoma povecujemo stevilo
       klastrov,
       dokler vsak klaster ne zadostuje dolocenim geom. kriterijem
    */
    if (znPP == 1)
        do_ligands(znPP, 1, PP_DIST, clusidPP, maskidPP, zdataPP, zpligPP, _pp);
    for (int k = 1; k < NCLTRIES && k < znPP; k++) {
        cuttree(znPP, NPP, k, clusidPP);
#ifdef VERB
        for (int i = 0; i < znPP; i++) {
            cout << "PP : it = " << k << " n = " << i
                 << " maskidPP = " << maskidPP[i] << " cl = " << clusidPP[i]
                 << endl;
        }
#endif
        do_ligands(znPP, k, PP_DIST, clusidPP, maskidPP, zdataPP, zpligPP, _pp);
    }
    //#ifdef VERB
    clock_t start2 = clock();
    //#endif
    if (znSL == 1)
        do_ligands(znSL, 1, SL_DIST, clusidSL, maskidSL, zdataSL, zpligSL, _sl);
    for (int k = 1; k < NCLTRIES && k < znSL; k++) {
        cuttree(znSL, NSL, k, clusidSL);
#ifdef VERB
        for (int i = 0; i < znSL; i++) {
            cout << "SL : it = " << k << " n = " << i
                 << " maskidSL = " << maskidSL[i] << " cl = " << clusidSL[i]
                 << endl;
        }
#endif
        do_ligands(znSL, k, SL_DIST, clusidSL, maskidSL, zdataSL, zpligSL, _sl);
    }

    if (znNU == 1)
        do_ligands(znNU, 1, NU_DIST, clusidNU, maskidNU, zdataNU, zpligNU, _nu);
    for (int k = 1; k < NCLTRIES && k < znNU; k++) {
        cuttree(znNU, NNU, k, clusidNU);
#ifdef VERB
        for (int i = 0; i < znNU; i++) {
            cout << "NU : it = " << k << " n = " << i
                 << " maskidNU = " << maskidNU[i] << " cl = " << clusidNU[i]
                 << endl;
        }
#endif
        do_ligands(znNU, k, NU_DIST, clusidNU, maskidNU, zdataNU, zpligNU, _nu);
    }

    if (znIO == 1)
        do_ligands(znIO, 1, IO_DIST, clusidIO, maskidIO, zdataIO, zpligIO, _io);
    for (int k = 1; k < NCLTRIES && k < znIO; k++) {
        cuttree(znIO, NIO, k, clusidIO);
#ifdef VERB
        for (int i = 0; i < znIO; i++) {
            cout << "IO : it = " << k << " n = " << i
                 << " maskidIO = " << maskidIO[i] << " cl = " << clusidIO[i]
                 << endl;
        }
#endif
        do_ligands(znIO, k, IO_DIST, clusidIO, maskidIO, zdataIO, zpligIO, _io);
    }

    //#ifdef VERB
    cout << "Ligands::cluster() Time of cuttree + do_ligands (SL) = "
         << (double)(clock() - start2) / CLOCKS_PER_SEC << endl;
    cout << "Ligands::cluster() Cumulative Time = "
         << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
    //#endif

    /* sprostimo pomnilnik in pocistimo :-)
     */
    for (int i = 0; i < nrows; i++) delete[] mask[i];
    delete[] mask;

    delete[] maskidPP;
    delete[] maskidSL;
    delete[] maskidNU;
    delete[] maskidIO;

    delete[] clusidPP;
    delete[] clusidSL;
    delete[] clusidNU;
    delete[] clusidIO;

    // C-stil, ker spomin allocatamo v clusterlib knjiznici
    free(NPP);
    free(NSL);
    free(NNU);
    free(NIO);
}

void Ligands::do_ligands(int nrows, int k, double radius, int* clusid,
                         int* maskid, double** data, Ligand** plig,
                         bstype _LIG) {
    /*
      Gremo po vseh klastrih, ki jih je vrgel ven cuttree in poiscemo tiste, ki
      ustrezajo (tocke v njih
      so dovolj blizu skupaj). Le te zapisemo kot bsite.
      Binding site score je izracunan po enacbi:  bs-score = (cluster-score)_max
      _max predtavlja ligand z najvisjim cluster-score-om.
     */

    /* gremo po vseh klastrih */
    for (int m = 0; m < k; m++) {
        bool do_check_cluster = false;
        // cluster (m) upostevamo samo, ce ni bil ze dodan kot binding site
        for (int i = 0; i < nrows; i++)
            if (clusid[i] == m && maskid[i] >= 0) {
                do_check_cluster = true;
                break;
            }

        /* ce je vec kot 90% tock manj kot PP_DIST (za protein-protein)
         * angstromov od sredisca klastra, je le ta ok!
         */
        Coor bsc;
        if (do_check_cluster && check_cluster(radius, clusid, m, nrows, data,
                                              bsc)) {  // v bsc je center bs-ja

            BSite* bs = new BSite(bsc, _LIG);
            bs->lig.clear();

#ifdef VERB
            cout << "Ligands::do_ligands()   binding site (x,y,z) = "
                 << bs->center.x << "," << bs->center.y << "," << bs->center.z
                 << endl;
#endif

            // oznacimo, da tega klastra (m-tega) ne delimo vec naprej in dodamo
            // ligande v bsite
            for (int i = 0; i < nrows; i++) {
                if (clusid[i] == m) {
                    maskid[i] = -1;
                    // dodamo ligande, ki pripadajo temu binding site-u
                    bs->lig.insert(plig[i]);
                    for (Ligand* li = plig[i]->first_neighb(); li != NULL;
                         li = plig[i]->next_neighb()) {
                        bs->lig.insert(li);
                    }
                }
            }

            /* dolocimo binding site score */
            // dobimo ligand z najvisjim cluster_score
            Ligand* limax =
                *max_element(bs->lig.begin(), bs->lig.end(), by_cluster_score);

            // izracunamo bs-score
            bs->score = limax->cd->cluster_score;

            // dodamo v unijo residueje na query proteinu, ki pripadajo temu
            // binding site-u (unija residuejev vseh ligandov)
            for (set<Ligand*>::iterator it = bs->lig.begin();
                 it != bs->lig.end(); it++) {
                // unikatno dodamo vse query residueje pripadajoce enemu ligandu
                // v unijo
                bs->unija.insert((*it)->rlist.begin(), (*it)->rlist.end());
            }
            // dodamo klaster kot en protein-protein binding site
            bsite.push_back(bs);
        }
    }
}

bool Ligands::check_cluster(double radius, int* clusid, int cl_id, int nrows,
                            double** data, Coor& center) {
    /*
       Izracunamo geometrijsko sredisce klastra tock in izracunamo, koliko jih
       pade v dolocen podan radij
       od sredisca, in koliko jih pade izven tega radija. Radij je razlicen za
       razlicne ligande, ki jih
       doloca mask. Ce je vec kot 90% tock znotraj radija, potem funkcija vrne
       true.
       Poleg tega je v center vrnjeno sredisce klastra
     */

    set_zero(center);

    int num_in_clus = 0;

    /* gremo po vseh pripadnikih klastra (cl_id) in izracunamo center klastra */
    for (int i = 0; i < nrows; i++) {
        if (clusid[i] == cl_id) {
            center.x += data[i][0];
            center.y += data[i][1];
            center.z += data[i][2];
            num_in_clus++;

#ifdef VERB
            cout << "binding site (x,y,z) = " << center.x << "," << center.y
                 << "," << center.z << endl;
#endif
        }
    }

    if (num_in_clus == 0) return false;  // je to mogoce ??

    /* izracunamo sredisce klastra */
    center = center * (1 / (double)num_in_clus);
    int num_good = 0;
    /* sedaj preverimo razdalje pripadnikov od centra */
    for (int i = 0; i < nrows; i++) {
        /* vzamemo samo tiste elemente, ki pripadajo cl_id-temu klastru */
        if (clusid[i] == cl_id) {
            Coor el;

            el.x = data[i][0];
            el.y = data[i][1];
            el.z = data[i][2];

            if (dist(center, el) < radius) num_good++;
        }
    }
#ifdef VERB
    cout << "num_good = " << num_good << " num_in_clus = " << num_in_clus
         << endl;
#endif
    return ((double)num_good / (double)num_in_clus > 0.9) ? true : false;
}

void Ligands::output() {
    /*
       Appendamo informacije o klastriranih ligandih v results.html - na konec.
     */

    //#ifdef VERB
    clock_t start = clock();
    //#endif

    //  ofstream out_file("ligands.json", ios_base::app);
    ofstream out_file((add_end_slash(OUTDIR) + "ligands.json").c_str());

    int bsnum = 0;

    //  sort(bsite.begin(), bsite.end(), by_bstype_bsscore);
    sort(bsite.begin(), bsite.end(), by_bstype_size);

    /* izpisemo se javascript objekt za binding site in ligande v JSON formatu
     */
    out_file << "[";

    for (vector<BSite*>::iterator it = bsite.begin(); it != bsite.end(); it++) {
        BSite* bs = *it;

        // vsak binding site z vsemi ligandi je en objekt v arrayu
        out_file << (it == bsite.begin() ? "{" : ",{");
        //    out_file << "\"num\":" << bsnum << ",\"type\":" << bs->type <<
        //    ",\"score\":" << bs->score << ",\"size\":" << size <<
        //    ",\"qresi\":[";
        out_file << "\"num\":" << bsnum << ",\"type\":" << bs->type
                 << ",\"score\":" << bs->score << ",\"qresi\":[";

        // izpisemo vse residueje, ki pripadajo dolocenemu binding site-u
        for (set<ChResi>::iterator it2 = bs->unija.begin();
             it2 != bs->unija.end(); it2++)
            out_file << ((it2 == bs->unija.begin()) ? "\"" : ",\"")
                     << it2->second << ":" << it2->first << "\"";

        // center binding site-a
        out_file << "],\"center\":{\"x\":" << bs->center.x
                 << ",\"y\":" << bs->center.y << ",\"z\":" << bs->center.z
                 << "},";

        // ligandi so objekti, ki pripadajo krovnemu objektu binding site-u
        out_file << "\"lig\":[";

        for (set<Ligand*>::iterator it2 = bs->lig.begin(); it2 != bs->lig.end();
             it2++) {
            Ligand* li = *it2;

            out_file << (it2 == bs->lig.begin() ? "{" : ",{");

            // izpisemo podatke za ligand ...
            //      out_file << "\"pdb_id\":\"" << li->pdb_id <<
            //      "\",\"cluster_id\":\"" << li->cd->cluster_id <<
            //      "\",\"chain_id\":\"" << li->chain_id
            out_file << "\"pdb_id\":\"" << li->pdb_id << "\",\"acid\":\""
                     << li->acid << "\",\"alignment_no\":\""
                     << li->cd->cluster_id << "\",\"chain_id\":\""
                     << li->chain_id  // NOVO
                     << "\",\"model\":\"" << li->model << "\",\"resn\":\""
                     << remove_whitespaces(li->resn)
                     << "\",\"resi\":" << li->resi
                     << ",\"score\":" << li->cd->cluster_score;

            // ... tudi ime rotirane ligand datoteke
            //      string fname = li->pdb_id + "." +
            //      to_string(li->cd->cluster_id) + "." + to_string(li->model) +
            //      "." + to_string(li->chain_id) + "."
            string fname = li->pdb_id + to_string(li->acid) + "." +
                           to_string(li->cd->cluster_id) + "." +
                           to_string(li->model) + "." +
                           to_string(li->chain_id) + "."  // NOVO
                           + remove_whitespaces(li->resn) + "." +
                           to_string(li->resi) + ".lig.pdb";

            out_file << ",\"f\":\"" << fname;

            // ... pa se kemijsko ime liganda oziroma molekule (ki ga moramo
            // prebrati iz ustrezne datoteke)
            string name = add_end_slash(LIGDIR) + "names/" +
                          ((bs->type == _pp || bs->type == _nu)
                               ? (li->pdb_id + li->chain_id)
                               : remove_whitespaces(li->resn));
            ifstream f(name.c_str());
            string text = "-";
            if (f.is_open()) {
                getline(f, text);
            }
            f.close();

            /* escapamo posebne znake : " in \ */
            text = json_escape_special_characters(text);

            out_file << "\",\"name\":\"" << text;

            out_file << "\",\"qresi\":[";
            // izpisemo resi-je in chain_id-je ligand binding site-ov
            for (set<ChResi>::iterator it3 = li->rlist.begin();
                 it3 != li->rlist.end(); it3++) {
                //        out_file << (( it3 == li->rlist.begin() ) ? "" : ",")
                //        << it3->second ;
                out_file << ((it3 == li->rlist.begin()) ? "" : ",") << "\""
                         << it3->second << ":" << it3->first << "\"";
            }

            //      // izpisemo chain_id-je ligand binding site-ov
            //      out_file << "],\"qcid\":[";
            //      for (set<ChResi>::iterator it3 = li->rlist.begin(); it3 !=
            //      li->rlist.end(); it3++) {
            //        out_file << (( it3 == li->rlist.begin() ) ? "\"" : ",\"")
            //        << it3->first << "\"";
            //      }

            out_file << "]}";
        }

        // zakljucimo lig array ligandov
        out_file << "]";
        // zakljucimo en binding site
        out_file << "}";
        bsnum++;
    }

    // zakljucimo array binding site-ov (to je ena __LIGS vrstica)
    out_file << "]" << endl;
    out_file.close();

    //#ifdef VERB
    cout << "Ligands::output() Time = "
         << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
//#endif

#ifdef VERB
    int i = 0;

    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        cout << "__BSIT:" << i << "#"
             << "IO"
             << "#" << 0 << "#" << 0 << "#" << (*it)->crd.x << "#"
             << (*it)->crd.y << "#" << (*it)->crd.z << endl;
        i++;
    }

    bsnum = 0;

    for (vector<BSite*>::iterator it = bsite.begin(); it != bsite.end(); it++) {
        BSite* bs = *it;

        string str;

        switch (bs->type) {
            case _pp:
                str = "PP";
                break;
            case _sl:
                str = "SL";
                break;
            case _nu:
                str = "NU";
                break;
            case _io:
                str = "IO";
                break;
        }

        for (set<Ligand*>::iterator it2 = bs->lig.begin(); it2 != bs->lig.end();
             it2++) {
            Ligand* li = *it2;

            // izpisemo kode ligandov za vsak binding site
            cout << "Ligands::output()  bsnum = " << bsnum
                 << " cluster_score = " << li->cd->cluster_score
                 << " pdb_id = " << li->pdb_id << " resn = "
                 << " acid = " << li->acid << " resn = " << li->resn
                 << " resi = " << li->resi
                 << " cluster_id = " << li->cd->cluster_id
                 << " x = " << li->crd.x << " y = " << li->crd.y
                 << " z = " << li->crd.z << endl;

            // izpisemo rotacijsko matriko in translacijski vektor
            cout << "Ligands::output() MATRIX "
                 << gsl_matrix_get(li->cd->U, 0, 0) << " "
                 << gsl_matrix_get(li->cd->U, 0, 1) << " "
                 << gsl_matrix_get(li->cd->U, 0, 2) << " "
                 << gsl_matrix_get(li->cd->U, 1, 0) << " "
                 << gsl_matrix_get(li->cd->U, 1, 1) << " "
                 << gsl_matrix_get(li->cd->U, 1, 2) << " "
                 << gsl_matrix_get(li->cd->U, 2, 0) << " "
                 << gsl_matrix_get(li->cd->U, 2, 1) << " "
                 << gsl_matrix_get(li->cd->U, 2, 2) << " " << endl;

            cout << "Ligands::output() VECTOR " << gsl_vector_get(li->cd->t, 0)
                 << " " << gsl_vector_get(li->cd->t, 1) << " "
                 << gsl_vector_get(li->cd->t, 2) << " " << endl;
        }
        bsnum++;
    }
#endif
}

void Ligands::write_rotate_ligands(Molecule* m, bool _caonly) {
    /*
       Izlocimo ligande iz .bu.pdb datotek (ki se nahajajo v $LIGDIR/biounits in
       vsakega zapisemo v svojo datoteko.
       Koordinate rotiramo, tako da so alignane na query. Tudi ce kaksna .bu.pdb
       datoteka manjka, ni problema - naredil bo
       za vse, ki obstajajo. Ce je vklopljen _caonly == true, bo izpisal samo CA
       atome v proteinskih ligandih.
    */

    //#ifdef VERB
    clock_t start = clock();
    //#endif

    Kabsch* k = new Kabsch();
    vector<Ligand*> l;
    /* ligande uredimo po pdb_id (oziroma pointerju na Molecule), nato po
     * cluster_id, nato po stevilki modela, chain_id in resi */
    for (set<Ligand*, ligcomp>::iterator it = lig_query.begin();
         it != lig_query.end(); it++) {
        l.push_back(*it);
    }

    sort(l.begin(), l.end(), by_ligcomp);  // NOVO

#ifdef VERB
    for (vector<Ligand*>::iterator it = l.begin(); it != l.end(); it++) {
        //    cout << "Ligands::write_ligands()  : sortedlig == " <<
        //    to_string((*it)->cd->m->pdb_id) + "." + (*it)->pdb_id + "." +
        //    to_string((*it)->cd->cluster_id) + "." + to_string((*it)->model) +
        //    "." + to_string((*it)->chain_id) + "." + (*it)->resn + "." +
        //    to_string((*it)->resi) + ".lig.pdb" << endl;
        cout << "Ligands::write_ligands()  : sortedlig == "
             << to_string((*it)->cd->m->pdb_id) + "." + (*it)->pdb_id +
                    to_string((*it)->acid) + "." +
                    to_string((*it)->cd->cluster_id) + "." +
                    to_string((*it)->model) + "." + to_string((*it)->chain_id) +
                    "." + (*it)->resn + "." + to_string((*it)->resi) +
                    ".lig.pdb"
             << endl;  // NOVO
    }
#endif

    int model = 0;
    string pdb_id = "";
    /* gremo po vseh sortiranih ligandih (a ne povecujemo it tukaj !!!!) */
    for (vector<Ligand*>::iterator it = l.begin(); it != l.end();) {
        /* dobimo vse ligande, ki pripadajo enakemu klastru (isti *Molecule in
         * cluster_id) */
        ClusterData* cd = (*it)->cd;
        vector<Ligand*>::iterator first = it;
        while (it != l.end() && (*it)->cd == cd) {
            it++;
        }
        vector<Ligand*>::iterator last = it;

#ifdef VERB
        cout << "Ligands::write_ligands()  : first == ";
        (*first)->output();
        if (*last) {
            cout << "Ligands::write_ligands()  : last == ";
            (*last)->output();
        } else {
            cout << "Ligands::write_ligands()  : last == l.end()" << endl;
        }
#endif

        /* odpremo ustrezno .bu.pdb datoteko (od predstojnika) */
        string rotated_name = add_end_slash(LIGDIR) + "biounits/" +
                              cd->m->pdb_id + cd->m->chain_id + ".bu.pdb";

        /* preverimo, ce ta datoteka sploh obstaja preden jo preberemo */
        ifstream frotated(rotated_name.c_str());
        if (!frotated.is_open()) {
            cout << "Warning (LIGANDS) : Missing " << rotated_name
                 << " ... Skipping!" << endl;
            continue;
        }

        /* odpremo datoteke v katere bomo izpisovali trenutne ligande */
        map<string, ofstream*> fligand;
        for (vector<Ligand*>::iterator it2 = first; it2 != last; it2++) {
            //      string ligand_name = (*it2)->pdb_id + "." +
            //      to_string((*it2)->cd->cluster_id) + "." +
            //      to_string((*it2)->model) + "." + to_string((*it2)->chain_id)
            //      + "." + (*it2)->resn + "." + to_string((*it2)->resi) +
            //      ".lig.pdb";
            //      string ligand_name = (*it2)->pdb_id + "." +
            //      to_string((*it2)->cd->cluster_id) + "." +
            //      to_string((*it2)->model) + "."
            string ligand_name = add_end_slash(OUTDIR) + (*it2)->pdb_id +
                                 to_string((*it2)->acid) + "." +
                                 to_string((*it2)->cd->cluster_id) + "." +
                                 to_string((*it2)->model) + "."  // NOVO
                                 + to_string((*it2)->chain_id) + "." +
                                 remove_whitespaces((*it2)->resn) + "." +
                                 to_string((*it2)->resi) + ".lig.pdb";
            ofstream* fl = new ofstream();
            fl->open(ligand_name.c_str());
            fligand[ligand_name] = fl;
        }

        char acid = ' ';  // NOVO

        /* ker so ligandi v l sortirani, se v datoteki nahajajo po vrsti, tako
         * da gremo le enkrat skozi */
        while (!frotated.eof()) {
            string s;
            getline(frotated, s);
            if (s.compare(0, 5, "MODEL") == 0) {
                model = atoi(s.substr(10, 4).c_str());
                //        pdb_id = s.substr(19,5);
                pdb_id = s.substr(19, 4);
                acid = s.at(23);  // NOVO
            }

            if (s.compare(0, 4, "ATOM") == 0 ||
                s.compare(0, 6, "HETATM") == 0) {
                /* ce naletimo na ligand */
                string resn = s.substr(17, 3);
                int resi = atoi(s.substr(22, 4).c_str());
                char chain_id = s.at(21);
                string tag = s.substr(13, 3);
                if (s.compare(0, 4, "ATOM") == 0 &&
                    amino.find(resn) !=
                        string::npos) {  // zagotoviti moras, da je ATOM
                                         // keyword, ker pod HETATM so lahko
                                         // tudi aminokisline (ki pa niso del
                                         // proteina)
                    resn = "@@@";
                    resi = 1;
                } else if (nucleic.find(resn) != string::npos) {
                    resn = "+++";
                    resi = 1;
                }

                // ce je proteinski ligand, potem izpisemo samo CA atome
                if (resn != "@@@" || tag == "CA " || !_caonly) {
                    //          Ligand *li = new Ligand(model, pdb_id, chain_id,
                    //          resn, resi, cd);
                    Ligand* li = new Ligand(model, pdb_id, acid, chain_id, resn,
                                            resi, cd);  // NOVO

#ifdef VERB
                    cout << "Ligands::write_ligands()  : li == ";
                    li->output();
#endif
                    if (binary_search(first, last, li, by_ligcomp)) {  // NOVO
                        //            string ligand_name = li->pdb_id + "." +
                        //            to_string(li->cd->cluster_id) + "." +
                        //            to_string(li->model) + "."
                        string ligand_name =
                            li->pdb_id + to_string(li->acid) + "." +
                            to_string(li->cd->cluster_id) + "." +
                            to_string(li->model) + "."  // NOVO
                            + to_string(li->chain_id) + "." +
                            remove_whitespaces(li->resn) + "." +
                            to_string(li->resi) + ".lig.pdb";
#ifdef VERB
                        cout << "Ligands::write_ligands()  : ligand_name = "
                             << ligand_name << endl;
#endif
                        /* rotiramo koordinate */
                        Coor c;

                        c.x = atof(s.substr(30, 8).c_str());
                        c.y = atof(s.substr(38, 8).c_str());
                        c.z = atof(s.substr(46, 8).c_str());
                        /* uporabiti moramo inverzno rotacijo, ker rotiramo
                         * alignan protein na query */
                        Coor crot =
                            k->rotate_vector_inv(c, li->cd->U, li->cd->t);
                        char temp[SMALL];

                        sprintf(temp, "%8.3f%8.3f%8.3f", crot.x, crot.y,
                                crot.z);
                        s.replace(30, 24, temp);

                        *fligand[ligand_name] << s << endl;
                    }

                    delete li;
                }
            }
        }  // while

        // zapremo ligand datoteke
        for (map<string, ofstream*>::iterator it2 = fligand.begin();
             it2 != fligand.end(); it2++) {
            it2->second->close();
        }

        frotated.close();
#ifdef VERB
        cout << "Ligands::write_ligands()  : rotated_name = " << rotated_name
             << endl;
#endif
    }

    delete k;

    //#ifdef VERB
    cout << "Ligands::write_rotate_ligands() Time = "
         << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
    //#endif
}

void Ligands::pdb_header(const string name) {
    /*
       Preberemo header sekcijo dane PDB datoteke, HETNAM in COMPDN kljucni
       besedi, in dobimo mapiranje med
       kodo liganda in njegovim kemijskim imenom ter pdb_idchain_id verige in
       njenim imenom.
    */
    ifstream f((add_end_slash(INDIR) + name).c_str());
    if (!f.is_open()) {
        cout << "Warning (LIGANDS) : Cannot open file for reading header "
             << name << " ... Skipping!" << endl;
        return;
    }

    /* preberemo samo header sekcijo v PDB datoteki */
    bool izstop = false;
    string text = "";
    bool berem_molekulo = false, berem_verigo = false;
    string chains = "";
    map<char, string> compnd;
    map<string, string> hetnam;

    while (!(f.eof() || izstop)) {
        string s;
        getline(f, s);
        /* ko pridemo do ATOM je konec */
        if (s.compare(0, 4, "ATOM") == 0 || s.compare(0, 6, "HETATM") == 0) {
            izstop = true;
        }

        /* opis za male ligande */
        if (s.compare(0, 6, "HETNAM") == 0) {
            /* ce ni t.i. continuation stevilke, potem smo na zapisu za nov
             * ligand */
            if (s.substr(9, 1) == " ") {
                hetnam[s.substr(11, 3)] = "";
            }
            // ... ce se vrstica ne konca z '-', potem dodamo presledek
            if (hetnam[s.substr(11, 3)].size() > 0) {
                if (hetnam[s.substr(11, 3)].at(hetnam[s.substr(11, 3)].size() -
                                               1) != '-')
                    hetnam[s.substr(11, 3)] += " ";
            }
            hetnam[s.substr(11, 3)] += remove_whitespaces(s.substr(15, 55));
        }

        /* opis za proteine in nukleinske kisline */
        if (s.compare(0, 6, "COMPND") == 0) {
            /* nov opis za molekulo */
            if (s.compare(11, 8, "MOLECULE") == 0) {
                text = remove_whitespaces(s.substr(21, 60));
#ifdef VERB
                cout << "Ligands::pdb_header()  : " << text << endl;
#endif
                berem_molekulo = true;
            } else if (s.compare(11, 5, "CHAIN") == 0) {
                chains = remove_whitespaces(s.substr(18, 63));
#ifdef VERB
                cout << "Ligands::pdb_header()  : " << chains << endl;
#endif
                berem_verigo = true;
            } else if (berem_molekulo) {
                if (text.size() > 0)
                    if (text.at(text.size() - 1) != '-')
                        text += " ";  // ... ce se vrstica ne konca z '-', potem
                                      // dodamo presledek (e.g. 1.vrsta: DNA
                                      // (5'- 2.vrsta: D(*TP*GP*T ...)
                text += remove_whitespaces(s.substr(11, 70));
#ifdef VERB
                cout << "Ligands::pdb_header()  : " << text << endl;
#endif
            } else if (berem_verigo) {
                chains += remove_whitespaces(s.substr(11, 70));
#ifdef VERB
                cout << "Ligands::pdb_header()  : " << chains << endl;
#endif
            }

            // ce najdemo kaksen drug keyword oz. podpicje v text ali chains,
            // potem zakljucimo branje
            if (text.find(';') != string::npos ||
                s.compare(11, 5, "CHAIN") == 0) {
                berem_molekulo = false;
            }
            if (chains.find(';') != string::npos ||
                s.compare(11, 8, "FRAGMENT") == 0 ||
                s.compare(11, 7, "SYNONYM") == 0 ||
                s.compare(11, 10, "ENGINEERED") == 0) {
                berem_verigo = false;
                // znebimo se nadleznega podpicja ;-)
                size_t found = text.find_last_of(';');
                if (found != string::npos) text.erase(found);
                // izlocimo iz chains posamezne chain_id-je (vejic in presledkov
                // ne smemo upostevati)
                for (string::iterator it = chains.begin(); it != chains.end();
                     it++) {
                    if (isalnum(*it)) {
                        compnd[*it] = text;
                    }
                }
                chains.clear();
            }
        }
    }
    cout << endl;

    // se zadnji chain_id ! (ce je na primer takoj za "COMPND    CHAIN: B" nek
    // drug keyword, e.g. "SOURCE ...." glej 1fss.pdb
    // znebimo se nadleznega podpicja ;-)
    size_t found = text.find_last_of(';');
    if (found != string::npos) text.erase(found);

    // izlocimo iz chains posamezne chain_id-je (vejic in presledkov ne smemo
    // upostevati)
    for (string::iterator it = chains.begin(); it != chains.end(); it++) {
        if (isalnum(*it)) {
            compnd[*it] = text;
        }
    }
    chains.clear();
    text.clear();

    /* izpisemo v datoteke */
    for (map<string, string>::iterator it = hetnam.begin(); it != hetnam.end();
         it++) {
        string out_name = add_end_slash(OUTDIR) + remove_whitespaces(it->first);
        ofstream fhet(out_name.c_str());
        fhet << it->second << endl;
        fhet.close();
#ifdef VERB
        cout << out_name << " --> " << it->second << endl;
#endif
    }
    string pdb_id = extract_pdb_code(name);
    pdb_id.erase(0, 3);  // zbrisemo 'pdb' iz imena 'pdb1all'
    for (map<char, string>::iterator it = compnd.begin(); it != compnd.end();
         it++) {
        string out_name = add_end_slash(OUTDIR) + pdb_id + it->first;
        ofstream fcom(out_name.c_str());
        fcom << it->second << endl;
        fcom.close();
        //#ifdef VERB
        cout << out_name << " --> " << it->second << endl;
        //#endif
    }

    f.close();
}
