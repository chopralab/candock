#include "molecule.h"
#include "desc.h"
#include "ligands.h"
#include "atom.h"
#include "probe.h"
#include "eelement.h"
#include "biounit.h"
#include "grid.h"
#include "kabsch.h"
#include "item.h"
#include "output.h"
#include "clusterdata.h"
#include "motif.h"
#include "debug.hpp"

/* callback funkcije */
bool by_flx_cluster_score(const pair<Residue*, Residue*> &i, const pair<Residue*, Residue*> &j) { 
  //  flx  cluster_score
  //  ------------------------------------
  //  true    0.1              true    0.1
  //  false   0.5     --->     true    2.0
  //  true    2.0              false   0.5
  //  false   1.5              false   1.5  ---> max_element
  //
  if (i.second->flx && !j.second->flx) 
    return true;
  else if ( i.second->flx == j.second->flx) 
    if (i.second->cd->cluster_score < j.second->cd->cluster_score) 
      return true;
  return false;
}

bool by_noflx_cluster_score(const pair<Residue*, Residue*> &i, const pair<Residue*, Residue*> &j) { 
  if (i.second->cd->cluster_score < j.second->cd->cluster_score) 
    return true;
  return false;
}



bool by_conservation(const pair<ChResi,int> &i, const pair<ChResi, int> &j) { return i.second > j.second; }

Molecule::Molecule(string n, string c_id, int m_id) {
#ifdef VERB
  cout << "Initializing molecule " << n << endl;
  //~ return;
#endif
  size_t found = n.find_last_of(".");  // najdemo suffix
  string suffix = n.substr(found + 1);
  if (!(suffix == "pdb"||suffix == "PDB"||suffix == "ent"||suffix == "ENT"||suffix == "srf"||suffix == "SRF")) {
    throw Err("Error (MOLECULE) : Input protein structure must be in Protein Data Bank format (.pdb or .ent) or surface format (.srf see --extract option)!", 10);
  }

  main_id = m_id;
  size = 0; 
  name = n;
  atom = NULL;
  comp = NULL;
  nfp = 0;
  pdb_id = extract_pdb_code(n);

  // ce web server uporablja probis, potem vzamemo samo prve stiri znake v pdb_id prileganega proteina
  if (!_longnames) {
    pdb_id = pdb_id.substr(0,4);
  }

  chain_id = c_id;

}

Molecule::~Molecule() {
#ifdef VERB
  cout << "Deleting molecule ..." << endl;
#endif
  Atom *tmp;
  while (atom != NULL) {
    tmp = atom;
    atom = atom->next;
    delete tmp;
  }
  while (comp != NULL) {
    tmp = comp;
    comp = comp->next;
    delete tmp;
  }
  
  /* izbrisemo listo prileganih molekul */
  for (vector<Molecule*>::iterator it=lista_molekul.begin(); it != lista_molekul.end(); it++) 
    delete *it;
  lista_molekul.clear();

  /* sprostimo pomnilnik za biounite - rotacijske matrike in vektorje deletamo v BioUnit destruktorju */
  for (map<pair<int,int>, BioUnit*>::iterator it = biounit.begin(); it != biounit.end(); it++) {
    delete it->second;
  }
  biounit.clear();

  /* sprostimo pomnilnik za align_score */  
  //  for (map<pair<int, int >, ClusterData* >::iterator it = align_score.begin(); it != align_score.end(); it++) {
  for (map<pair<int, int>, ClusterData* >::iterator it = align_score.begin(); it != align_score.end(); it++) {
    delete it->second;
  }
  align_score.clear();

  /* sprostimo pomnilnik za sequence */  
  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
    delete *it;
  }
  sequence.clear();

  /* sprostimo pomnilnik za align tabelo, v kateri je povezava med prilegano ak in query ak */  
  for (multimap<Residue*, Residue*>::iterator it = align.begin(); it != align.end(); it++) {
    delete it->second;  // sprostimo pomnilnik za ak v alignanem proteinu (v query proteinu smo ze - glej sequence)
  }
  align.clear();
//  /* sprostimo pomnilnik za plain_pdb */  
//  for (vector<string*>::iterator it = plain_pdb.begin(); it != plain_pdb.end(); it++) {
//    if (*it) delete *it;  // sprostimo pomnilnik za ak v alignanem proteinu (v query proteinu smo ze - glej sequence)
//  }
}
 

bool Molecule::check_next_resi(size_t i) { 
  /*
    Preverimo, ali je resi naslednje aminokisline manjsi od trenutnega. Ce je, potem vrnemo false.
  */
  size_t j = i;
  while (j < plain_pdb.size() && plain_pdb[j].at(21) == plain_pdb[i].at(21) && atoi(plain_pdb[j].substr(22,4).c_str()) == atoi(plain_pdb[i].substr(22,4).c_str())) j++;
  if (j < plain_pdb.size() && plain_pdb[j].compare(0, 4, "ATOM") == 0 && plain_pdb[j].at(21) == plain_pdb[i].at(21) && atoi(plain_pdb[j].substr(22,4).c_str()) < atoi(plain_pdb[i].substr(22,4).c_str())) // JUL/06/2011 moras preveriti za ATOM, ker lahko je HETATM z nizjo stevilko (primer 1xyzB, MET 835)
    return false;
  return true;
}


void Molecule::read_PDB(bool _firstModel, int atomType, bool _allChains, int saveMem) {
  /* 
     Format PDB zapisa uposteva standarde na http://www.wwpdb.org/docs.html
  */
//  cout << "Molecule::read_PDB  :   At The Start!" << endl;
  
  int model = 1;
  int resi = -XLARGE;   // residue numbers in some pdbs can be < 0
  char tmp_chain_id = ' ';
  Atom *a = NULL, *la = NULL, *acid_start = NULL;
  int aasize = 0;
  int bio_num = 0; // stevilka bioloskega unita
//  int bio_cid = -1; // glej na primer PDB ID 1fss, kjer imas pod eno BIOMOLECULE vec APPLY THE FOLLOWING TO CHAINS
  set<char> bio_cid;

  /* gremo po plain PDB-ju, ki ga moramo prej prebrati z read_plain_pdb */
  for (size_t i = 0; i < plain_pdb.size(); i++) {
    

    /* preberemo samo prvi model, ali pa vse modele, pri izracunu interface residujev od prileganih ligandov*/
    if (_firstModel && plain_pdb[i].compare(0, 3, "END") == 0) {
      break;
    }
    
    /* preberemo zapis za vsak bioloski unit posebej */
    else if (saveMem != _saveMemLig && plain_pdb[i].compare(0, 10, "REMARK 350") == 0) {
      if (plain_pdb[i].compare(0, 23, "REMARK 350 BIOMOLECULE:") == 0) {
        bio_num++;
//        bio_cid = -1;
        pair <int, int> bio_id = make_pair(model, bio_num);
        biounit[ bio_id ] = new BioUnit();
      }
      
      /* najprej preberemo chain id-je, na katere se transformacije nanasajo */
      int reset = plain_pdb[i].compare(0, 41, "REMARK 350 APPLY THE FOLLOWING TO CHAINS:");
      if (reset == 0 ||
//      if (plain_pdb[i].compare(0, 41, "REMARK 350 APPLY THE FOLLOWING TO CHAINS:") == 0 ||
//          plain_pdb[i].compare(0, 41, "REMARK 350                    AND CHAINS:") == 0) {
          plain_pdb[i].find("AND CHAINS:") != string::npos) {
        if (reset == 0) bio_cid.clear();
        size_t found;
	
        found=plain_pdb[i].find_first_of(",");
        
        while (found!=string::npos) {
          bio_cid.insert( plain_pdb[i].at(found-1) ); 
          found=plain_pdb[i].find_first_of(",",found+1);
        }
        
        /* se zadnji chain id */
        size_t lci = plain_pdb[i].find_last_not_of(" \t\f\v\n\r");
        
        if ( isalnum(plain_pdb[i].at (lci) ) ) {
          bio_cid.insert( plain_pdb[i].at(lci) );
        }
      }
      
      /* preberemo rotacijske matrike U in translacijske vektorje t za ta biounit*/
      /* REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 */
      if (plain_pdb[i].compare(0, 18, "REMARK 350   BIOMT") == 0) {
        int biomt_num = atoi(plain_pdb[i].substr(18,1).c_str());
        int rota_num =  atoi(plain_pdb[i].substr(19,4).c_str());
        pair <int, int> bio_id = make_pair(model, bio_num);
	
	if (biounit[bio_id]) { // pri nekaterih pdb-jih (1oln) manjka REMARK 350 BIOMOLECULE, zato je prislo prej do seg. faulta
	  /* alociramo pomnilnik za U in t */
	  if (biomt_num == 1) {
	    biounit[bio_id]->U[rota_num] = gsl_matrix_alloc(3,3);
	    biounit[bio_id]->t[rota_num] = gsl_vector_alloc(3);
	    
	    // chain_id je vezan na rota_num
	    for (set<char>::iterator it = bio_cid.begin(); it != bio_cid.end(); it++) {
	      biounit[bio_id]->chain_id.insert(pair<int,char>(rota_num, *it));
	    }
	  }
	  
	  gsl_matrix_set(biounit[bio_id]->U[rota_num], biomt_num - 1, 0, atof(plain_pdb[i].substr(23,10).c_str()) );
	  gsl_matrix_set(biounit[bio_id]->U[rota_num], biomt_num - 1, 1, atof(plain_pdb[i].substr(33,10).c_str()) );
	  gsl_matrix_set(biounit[bio_id]->U[rota_num], biomt_num - 1, 2, atof(plain_pdb[i].substr(43,10).c_str()) );
	  
	  gsl_vector_set(biounit[bio_id]->t[rota_num], biomt_num - 1, atof(plain_pdb[i].substr(53,15).c_str()) );
	}
      }
    }
    else if (plain_pdb[i].compare(0, 6, "ENDMDL") == 0) {
      
      /* za vsak nov model resetiramo stevilko biomolekule */
      bio_num = 0;
              
      /* ce nimamo REMARK 350 zapisov, recimo, prebrali smo NMR strukturo, generiramo en biounit z 
         identiteto za matriko */
      pair <int, int> bio_id = make_pair(model, 1);
      
      if (biounit.find( bio_id ) == biounit.end() ) {
        
        biounit[ bio_id ] = new BioUnit();
        
        /* ne dodamo nobenega chaina v biounit[ bio_id ]->chain - to je znak, da se matrika nanasa na vse chaine */
    
        /* alociramo pomnilnik za U in t */
        biounit[ bio_id ]->U[1] = gsl_matrix_alloc(3,3);
        biounit[ bio_id ]->t[1] = gsl_vector_alloc(3);
    
        gsl_matrix_set_identity( biounit[ bio_id ]->U[1] );
        gsl_vector_set_zero( biounit[ bio_id ]->t[1] );
      }

    }

    else if (!_firstModel && plain_pdb[i].compare(0, 5, "MODEL") == 0) {

      /*      
              COLUMNS        DATA  TYPE    FIELD          DEFINITION
              ---------------------------------------------------------------------------------------
              1 -  6        Record name   "MODEL "
              11 - 14        Integer       serial         Model serial number.
              20 - 24        String        ime            PDB ID CHAIN ID (e.g. 1gotA)
              
              012345678901234567890123
              MODEL        1     1gotA
      */


      /* preberemo malo spremenjen MODEL zapis, vsak model (zacnejo se z 1) je en prilegan PDB file */
      model = atoi(plain_pdb[i].substr(10,4).c_str());

      modelName[model] = plain_pdb[i].substr(19,5);


    }
    /* preberemo MODRES zapise, ce obstajajo */
    else if (plain_pdb[i].compare(0, 6, "MODRES") == 0) {
      /*
        COLUMNS        DATA TYPE     FIELD       DEFINITION
        --------------------------------------------------------------------------------
        1 -  6        Record name   "MODRES"
        8 - 11        IDcode        idCode      ID code of this entry.
        13 - 15        Residue name  resName     Residue name used in this entry.
        17             Character     chainID     Chain identifier.
        19 - 22        Integer       seqNum      Sequence number.
        23             AChar         iCode       Insertion code.
        25 - 27        Residue name  stdRes      Standard residue name.
        30 - 70        String        comment     Description of the residue modification.
      */

      /* kot ligande ne upostevamo tistih MODRES-ov, ki so spremenjene normalne aminokisline, npr. MET --> MSE, in glikoziliranih aminokislin */

      /* samo, ce je insertion code blank */
      if (plain_pdb[i].at(22) == ' ')  {
        
        
        string stdRes = plain_pdb[i].substr(24, 3);
        
        /* naredimo seznam tistih spremenjenih ak, ki jih nocemo: ce je standardna residue v MODRES zapisu enaka eni izmed standardnih aminokislin */
        if ( amino.find(stdRes) != string::npos ) {
          
          
          int modResi = atoi( plain_pdb[i].substr(18, 4).c_str() );
          char ch     = plain_pdb[i].at(16);

#ifdef VERB
          cout << "read_PDB() : nasel modified residue s standardnim imenom = " << modResi << " " << ch << endl; 
#endif
          /* shranimo resi in chain id od spremenjenega residue-ja */
//          modres.insert( Ligand(model, "", ch, "", modResi, NULL) );
          modres.insert( Ligand(model, ch, modResi) );

        }
      }

    }

    else if (plain_pdb[i].compare(0, 6, "SITE  ") == 0) {
      /*
        COLUMNS        DATA  TYPE    FIELD         DEFINITION
        ---------------------------------------------------------------------------------
        1 -  6        Record name   "SITE  "
        8 - 10        Integer       seqNum        Sequence number.
        12 - 14        LString(3)    siteID        Site name.
        16 - 17        Integer       numRes        Number of residues that compose the site.
        19 - 21        Residue name  resName1      Residue name for first residue that 
        creates the site.
        23             Character     chainID1      Chain identifier for first residue of site.
        24 - 27        Integer       seq1          Residue sequence number for first residue
        of the  site.
        28             AChar         iCode1        Insertion code for first residue of the site.
        30 - 32        Residue name  resName2      Residue name for second residue that 
        creates the site.
        34             Character     chainID2      Chain identifier for second residue of
        the  site.
        35 - 38        Integer       seq2          Residue sequence number for second
        residue of the site.
        39             AChar         iCode2        Insertion code for second residue
        of the  site.
        41 - 43        Residue name  resName3      Residue name for third residue that 
        creates  the site.
        45             Character     chainID3      Chain identifier for third residue
        of the site.
        46 - 49        Integer       seq3          Residue sequence number for third
        residue of the site.
        50             AChar         iCode3        Insertion code for third residue
        of the site.
        52 - 54        Residue name  resName4      Residue name for fourth residue that 
        creates  the site.
        56             Character     chainID4      Chain identifier for fourth residue
        of the site.
        57 - 60        Integer       seq4          Residue sequence number for fourth
        residue of the site.
        61             AChar         iCode4        Insertion code for fourth residue
        of the site.

                 1         2         3         4         5         6         7         8
        12345678901234567890123456789012345678901234567890123456789012345678901234567890
        SITE     1 AC1  3 HIS A  94  HIS A  96  HIS A 119                               
        SITE     1 AC2  5 ASN A  62  GLY A  63  HIS A  64  HOH A 328                    
        SITE     2 AC2  5 HOH A 634                                                     
        SITE     1 AC3  5 GLN A 136  GLN A 137  PRO A 138  GLU A 205                    
        SITE     2 AC3  5 CYS A 206                                                     
        SITE     1 AC4 11 HIS A  64  HIS A  94  HIS A  96  HIS A 119                    
        SITE     2 AC4 11 LEU A 198  THR A 199  THR A 200  TRP A 209                    
        SITE     3 AC4 11 HOH A 572  HOH A 582  HOH A 635           

        Za opis glej REMARK800 zapis!        
      */

      string siteName = plain_pdb[i].substr(11, 3);

      int offRI = 23, offC = 22, offI = 27, offRN = 18;

      int siteResi = atoi ( plain_pdb[i].substr(offRI, 4).c_str() );
      string siteResn = plain_pdb[i].substr(offRN, 3);
      char ch = plain_pdb[i].at(offC);
      
      while (ch != ' ') {
        /* dodamo samo, ce je insertion code blank in je standardna aminokislina */
        if (plain_pdb[i].at(offI) == ' ' && amino.find(siteResn) != string::npos) 
          site.insert( make_pair( Ligand(model, ch, siteResi), siteName ) );

        offRN += 11;
        offRI += 11;
        offC += 11;
        offI += 11;

        siteResi = atoi ( plain_pdb[i].substr(offRI, 4).c_str() );
        siteResn = plain_pdb[i].substr(offRN, 3);
        ch = plain_pdb[i].at(offC);

      }



    }
    /* preberemo HETATM oziroma ATOM zapis */
    else if (
              ((atomType == _hetero || atomType == _nucleic) && plain_pdb[i].compare(0, 6, "HETATM") == 0 &&
               /* modified standardnih aminokislin ne upostevamo */
               modres.find( Ligand(model, plain_pdb[i].at(21), atoi( plain_pdb[i].substr(22, 4).c_str() ) ) ) == modres.end()
               
               ) 
              
              || 

              ( (plain_pdb[i].compare(0, 6, "HETATM") == 0 || plain_pdb[i].compare(0, 4, "ATOM") == 0) && 
                atomType == _nucleic && nucleic.find(plain_pdb[i].substr(17,3)) != string::npos 
              ) 

              ||

              ( plain_pdb[i].compare(0, 4, "ATOM") == 0 && 

                  /* residue name je ena izmed standardnih aminokislin */
                  amino.find(plain_pdb[i].substr(17,3)) != string::npos &&

                  /* tudi tu preverimo, da ni to modified standardna aminokislina, ki ima isto ime kot normalna, npr. TYR -> TYR */
                  modres.find( Ligand(model, plain_pdb[i].at(21), atoi( plain_pdb[i].substr(22, 4).c_str() ) ) ) == modres.end() &&
                  
                  /* preberemo samo query verigo oziroma vse verige */
                  ( _allChains || chain_id.find(plain_pdb[i].at(21)) != string::npos ) &&
                  /* zahtevamo, da je prvi znak atom tag-a prazen (kalcij "CA  " bo tako prebran kot "A  ", vodikov ne maramo */
                  plain_pdb[i].at(12) == ' ' && plain_pdb[i].at(13) != 'H' &&
                  
                  /* alternate location indicator */             
                  ( plain_pdb[i].at(16) == ' ' || plain_pdb[i].at(16) == 'A' ) &&

                  /* insertion code mora biti prazen - to so ak, ki so umetno dodane v verigo - da ne pokvarijo zaporedja, dodajo tu A */
                  plain_pdb[i].at(26) == ' ' &&
                  
                  /* residue number ne sme kar poskocit iz 188 na 1188 (in potem spet pasti na 189) - primer 1k9oE */
                  ( a == NULL || check_next_resi(i) )
              )
             ) {

//      if (saveMem != _saveMemBiou) {
      
      a = new Atom(size);
      
      if (comp == NULL) { comp = a; la = a; } else { la->next = a; la = a; }
      
      /* preberemo ATOM ali HETATM */
      // prvi znak (12) izpustimo, zato kalcij, "CA  ", pride kot "A  "
      
      a->tag[3] = '\0';
      
      strcpy(a->tag, plain_pdb[i].substr(13, 3).c_str() );
      
      
      if (a->tag[1] == ' ') a->tag[1] = '\0';
      if (a->tag[2] == ' ') a->tag[2] = '\0';
      
      
      strcpy(a->resn, plain_pdb[i].substr(17, 3).c_str() );
      
      a->chain_id = plain_pdb[i].at(21);
      a->resi = atoi( plain_pdb[i].substr(22, 4).c_str() );
      a->crd.x = atof( plain_pdb[i].substr(30, 8).c_str() );
      a->crd.y = atof( plain_pdb[i].substr(38, 8).c_str() );
      a->crd.z = atof( plain_pdb[i].substr(46, 8).c_str() );
      
      /* ce je HETATM, potem oznacimo v hetero */
      if (plain_pdb[i].compare(0, 6, "HETATM") == 0) a->het = true; else a->het = false;
      
      /* stevilka modela */
      a->model = model;
      
      
#ifdef VERB
      cout << "Molecule read_PDB " << a->tag << " " << a->resn << " " << a->chain_id << " " << a->resi 
           << " " << a->crd.x << " " << a->crd.y << " " << a->crd.z << " " << a->het << " " << a->model << endl;
#endif
      
      if (a->resi != resi || a->chain_id != tmp_chain_id) {
        tmp_chain_id = a->chain_id;
        resi = a->resi;
        if (acid_start != NULL) { 
          acid_start->aasize = aasize; 
          /* velikost aminokisline preverimo le, ce ni HETATM */
          if (!acid_start->het) acid_start->check_aasize(); 
          
          aasize = 0; 
        }
        
        acid_start = a;
      }
      
      a->acid_start = acid_start;
      
      switch(a->tag[0]) {
        
      case 'C': a->r = CATOM;    break; 
      case 'N': a->r = NATOM;    break;
      case 'O': a->r = OATOM;    break;
      case 'S': a->r = SATOM;    break;
      case 'H': a->r = HATOM;    break;
      default : a->r = CATOM;    break;     // ce atom ne pase v nobeno od kategorij, dobi radij C atoma
      }
      
      a->num = size;
      
      size++;
      aasize++;
//    } // konec odseka if _saveMemBiou ...
    }
  }  // konec plain_pdb zanke
  
  /* se zadnja prebrana aminokislina */
  if (acid_start != NULL) { 
    
    acid_start->aasize = aasize; 
    
    /* velikost aminokisline preverimo le, ce ni HETATM */
    if (!acid_start->het) acid_start->check_aasize(); 
    
    aasize = 0; 
    
  }
  
  int i = 0;
  int non_hetero = 0;
  
  /* preverimo, ali ni ves protein (>90% atomov pod ATOM) v cudnih aminokislinah */
  for (Atom *atm = comp; atm != NULL; atm = atm->next) {
    
    if (atm->acid_start->aasize == -1) i++;
    
    /* hetero atomov ne upostevamo, ampak oni tudi aasize nikoli nimajo -1 !!, tako da dve vrstici zgoraj HET ni treba preverjati */
    if (!atm->het) non_hetero++;
    
  }
  
  if ( i == non_hetero || i > (0.9 * non_hetero) ) {
    

    throw Err("Error (MOLECULE) : File " + name + " does not contain a protein: " + to_string(i) 
	      + " defect atoms were found out of " + to_string(non_hetero) + ".", 2);
  }
  
  /* preverimo, ce so (samo!) v prvem modelu kaksne prekrivajoce se koordinate */
  check_same_coor(comp);
  
#ifdef VERB
  for (map<pair<int,int>, BioUnit*>::iterator it = biounit.begin(); it != biounit.end(); it++) {
    
    int model = it->first.first;
    int bio_num = it->first.second;
    
    pair <int, int> bio_id;
    
    bio_id = make_pair(model, bio_num);

    cout << "read_PDB() MODEL " << model << " BIOMOLECULE " << bio_num << endl;
    for (set<char>::iterator it2 = it->second->chain.begin(); it2 != it->second->chain.end(); it2++) {
      cout << "read_PDB() CHAIN " << *it2 << endl;
    }

    for (map<int, gsl_matrix*>::iterator it2 = it->second->U.begin(); it2 != it->second->U.end(); it2++) {
      cout << "read_PDB() MATRIX " << it2->first << endl;
      cout << "read_PDB() MATRIX " << gsl_matrix_get(it2->second, 0, 0) << " " 
           << gsl_matrix_get(it2->second, 0, 1) << " " << gsl_matrix_get(it2->second, 0, 2) << " " << endl;
      cout << "read_PDB() MATRIX " << gsl_matrix_get(it2->second, 1, 0) << " " 
           << gsl_matrix_get(it2->second, 1, 1) << " " << gsl_matrix_get(it2->second, 1, 2) << " " << endl;
      cout << "read_PDB() MATRIX " << gsl_matrix_get(it2->second, 2, 0) << " " 
           << gsl_matrix_get(it2->second, 2, 1) << " " << gsl_matrix_get(it2->second, 2, 2) << " " << endl;
    }
    for (map<int, gsl_vector*>::iterator it2 = it->second->t.begin(); it2 != it->second->t.end(); it2++) {
      cout << "read_PDB() VECTOR " << gsl_vector_get(it2->second, 0) << " " 
           << gsl_vector_get(it2->second, 1) << " " << gsl_vector_get(it2->second, 2) << " " << endl;
    }
    
  }
#endif

}


void Molecule::check_same_coor(Atom *start) {
  /* 
     Preveri, od zacetka proteina, pa do konca, da katera dva atoma (v isti aminokislini) nimata istih koordinat.
     Ce najdemo enake koordiante, zapisemo v acid_start->aasize = -1, kar pomeni, da z aminokislino nekaj
     ni ok!

     UPDATE MAY/1/2011: Preverimo za dve zaporedni aminokislini (primer PDB ID 1oxk, chain A in C)
     UPDATE DEC/6/2009: To naredimo samo za ne hetero atome modela 1 !
  */

  for (Atom *at = start; at != NULL; at = at->next) {
//    for (Atom *at2 = at->next; at2 != NULL && at2->acid_start == at->acid_start; at2 = at2->next) {
    for (Atom *at2 = at->next; at2 != NULL && at2->resi <= at->resi + 1; at2 = at2->next) {
      if (!at->het && !at2->het && at->model == 1 && at2->model == 1 && dist_fast(at->crd, at2->crd) < EPS) {
        at->acid_start->aasize = -1;
        cout << "Warning (MOLECULE) : Residue " << at->resn << " " << at->resi << " contains atoms with identical coordinates!" << endl;
        
      }
    }
  }
}


void Molecule::restore_atoms() {
  /*
    For database search we have stored the surfaces (atoms and descriptors) of proteins. Here we read atoms
    directly into molecule->atom array (not the comp).
  */
  ifstream surf ((add_end_slash(INDIR) + name).c_str());
  char buffer[XLARGE];
  float r, x, y, z;
  int resi = -XLARGE;
  char tmp_chain_id = ' ';
  Atom *a, *la = NULL, *acid_start = NULL;
  
  if (!surf.is_open()) { 
    throw Err("Error (MOLECULE) : Cannot open SRF file " + name + ".", 1);
  }
  size = 0;
  while (!surf.eof() ) {
//    string str_buff;
//    getline(surf, str_buff);
//    char *buffer = new char[str_buff.size()];
    surf.getline (buffer,XLARGE);
//    cout << buffer << endl;
    if (strncmp(buffer, "A>", 2) == 0) {
      a = new Atom(0);
      if (atom == NULL) { atom = a; la = a; } else { la->next = a; la = a; }
//      sscanf(&buffer[2], "%5d%5f%8f%8f%8f%3d%6s%3s%5d", &a->num, &r, &x, &y, &z, &a->helix, a->tag, a->resn, &a->resi);

      /* helix je OBSOLETE, namesto njega v SRF file shranimo stevilko modela, ki jo zaenkrat ignoriramo (=1) */
//      sscanf(&buffer[2], "%5d%5f%8f%8f%8f%3d%6s%3s%5d", &a->num, &r, &x, &y, &z, &a->model, a->tag, a->resn, &a->resi);
      sscanf(&buffer[2], "%5d%5f%8f%8f%8f%*3d%6s%3s%5d", &a->num, &r, &x, &y, &z, a->tag, a->resn, &a->resi);

      a->r = r;
      a->crd.x = x;
      a->crd.y = y;
      a->crd.z = z;
      a->chain_id = buffer[54];

      a->model = 1; // stevilko modela postavimo na 1 zaradi minmax_coor() f-je (sicer jo tudi ze v Atom konstruktorju)


//      strncpy(a->chain_id, &buffer[53], 1);
//      a->chain_id[1]='\0';
      sscanf(&buffer[54], "%3d", &a->colors.size);
      for (int i=0; i < a->colors.size; i++)
        sscanf(&buffer[i*4 + 57], "%5d", &a->colors.color[i]);
      if (a->resi != resi || a->chain_id != tmp_chain_id) {
        resi = a->resi;
        tmp_chain_id = a->chain_id;
        acid_start = a;
      }
      a->acid_start = acid_start;
      //      a->print_atom();
      size++;
    }
    else if (strncmp(buffer, "MOL>", 4) == 0) {
      sscanf(&buffer[4], "%6d", &lcolor);
    }
  }
  surf.close();
  if (!size)
    throw Err("Error (MOLECULE) : Defect SRF file " + name + ".", 101);
}

void Molecule::restore_descriptors(Descriptor *&desc, bool _nobb) {
  /*
    For database search we have stored the surfaces (atoms and descriptors) of proteins. Here we read descriptors 
    into desc2 array. We have to set the desc->atom to the pointer to the correct atom (according to the read value of 
    atom->num).
    
    Ce je vklopljen _nobb, potem ne upostevamo deskriptorjev, ki kazejo na enega izmed backbone atomov: O, N
  */
  ifstream surf ((add_end_slash(INDIR) + name).c_str());
  char buffer[XLARGE];
//  float r, x, y, z;

  char c_id;

  Descriptor *d;
  Atom *a;
  desc = NULL;
  if (!surf.is_open()) { 
    throw Err("Error (MOLECULE) : Cannot open SRF file " + name + ".", 1);
  }
  //~ cout << "opening file for descriptors = " << name << endl;
  while (!surf.eof() ) {
    surf.getline (buffer, XLARGE);
    if (strncmp(buffer, "D>", 2) == 0) {
      d = new Descriptor(0);
      //      if (desc == NULL) { desc = d; ld = d; } else { ld->next = d; ld = d; }
      sscanf(&buffer[2], "%5d%5lf%8lf%8lf%8lf%3d%1d%5d%*2c%c", &d->num, &d->r, &d->crd.x, &d->crd.y, &d->crd.z, &d->s->mnsp, &d->psurf, &d->atom_num, &c_id);


      /* Inicializiramo d->atom */
      a = atom;
      
      // najdemo atom, ki pripada temu deskriptorju
      while (!(a->num == d->atom_num && a->chain_id == c_id)) a = a->next;
      /* ne upostevamo backbone deskriptorjev */
      if (_nobb) {
        /* ce je backbone, izbrisemo deskriptor in preskocimo nadaljne vrstice */
        if(strcmp(a->tag, "N") == 0 || strcmp(a->tag, "O") == 0) { //PEP 
          delete d;
          continue;
        }
      }
      
      d->atom = a;
      //    d->print_desc();
      
      /* dodamo deskriptor */
      d->next = desc;
      desc = d;
    }
  }
  surf.close();
  if (!desc)
    throw Err("Error (MOLECULE) : Defect SRF file " + name + ".", 101);
}


void Molecule::restore_probes(Probe *&probe) {
  /*
    For database search we have stored the surfaces (atoms, descriptors, and probes) of proteins. Here we read 
    probes into probe2 array.
  */
  ifstream surf ((add_end_slash(INDIR) + name).c_str());
  char buffer[XLARGE];
  float r, x, y, z;
  Probe *p;
  probe = NULL;

  if (!surf.is_open()) { 
    throw Err("Error (MOLECULE) : Cannot open SRF file " + name + ".", 1);
  }
  while (!surf.eof() ) {
    surf.getline (buffer, XLARGE);
    if (strncmp(buffer, "P>", 2) == 0) {
      p = new Probe(0);
      //      if (desc == NULL) { desc = d; ld = d; } else { ld->next = d; ld = d; }
      sscanf(&buffer[2], "%5d%5f%8f%8f%8f%5d", &p->num, &r, &x, &y, &z, &p->color);
      p->r = r;
      p->crd.x = x;
      p->crd.y = y;
      p->crd.z = z;
      p->next = probe;
      probe = p;
    }
  }
  surf.close();
  if (!probe)
    throw Err("Error (MOLECULE) : Defect SRF file " + name + ".", 101);
}


void Molecule::append_coor(string chain_id) {
  /* 
     We append (or first time initialize) the atoms from the comp to the atom array.
     Output is atom[].crd, atom[].r, atom[].resn, atom[].acid_start;
  */
  Atom *a, *la = NULL;
  Atom *tmp, *tmp1;
  /* First time initialization */
  if (atom == NULL) {
    tmp1 = comp;
    size = 0;
    while (tmp1 != NULL) {
      if ((chain_id.find(tmp1->chain_id) != string::npos || chain_id.empty()) && tmp1->acid_start->aasize > 0) {

        a = new Atom(size);
        if (atom == NULL) { atom = a; la = a; } else { la->next = a; la = a; }
        *a = *tmp1; // uporaba posebnega operatorja (glej atom.h)
        size++;
      }
      tmp1 = tmp1->next;
    }
  }
  /* Appendamo koordinate, uporabljamo samo pri state1 (PPI) */
  else {

#ifdef VERB
    cout << "Molecule::append_coor Appendamo koordinate !" << endl; 
#endif

    tmp = atom;
    while (tmp->next != NULL) tmp = tmp->next;
    tmp1 = comp;
    while (tmp1 != NULL) {
      if ((chain_id.find(tmp1->chain_id) != string::npos || chain_id.empty()) && tmp1->acid_start->aasize > 0) {
        a = new Atom(size);
        *a = *tmp1; // uporaba posebnega operatorja (glej atom.h)
        size++;
        tmp->next = a;
        tmp = a;
      }
      tmp1 = tmp1->next;
    }
  }
}


void Molecule::append_coor_delete_comp(string chain_id) {
  /* 
     Prepisemo vse atome iz comp v atom in hkrati brisemo comp, zato da porabimo manj pomnilnika pri
     ligandih, ko je vklopljena opcija --lig.
  */
  Atom *a, *la = NULL;
  Atom *tmp1;
  /* ne moremo appendat, samo ce je atom prazen */
  if (atom == NULL) {
    tmp1 = comp;
    size = 0;
    while (tmp1 != NULL) {
      if ((chain_id.find(tmp1->chain_id) != string::npos || chain_id.empty()) && tmp1->acid_start->aasize > 0) {

        a = new Atom(size);
        if (atom == NULL) { atom = a; la = a; } else { la->next = a; la = a; }
        *a = *tmp1; // uporaba posebnega operatorja (glej atom.h)
        size++;
      }
      Atom *prev = tmp1;
      tmp1 = tmp1->next;
      delete prev;
    }
  }

  comp = NULL;
}


void Molecule::init_grid_molecule(Grid *grid, float distance) {
  /* 
     Here the atomic coordinates are mapped to box-shaped segments. First, the atom with 
     the minimum x-coordinate is calculated. Second, the atoms are mapped to grid according
     to their relative coordinates.
     Note: First we delete the neighbor list for each atom, which is neccessary for 
     reinitializing the neighbor lists in state = 1.

     NEW (DEC/18/2009) : preden klices to funkcijo, mora biti inicializirana min in max
     koordinata (glej f-jo minmax_coor)
     Input is atom[].crd;
     Output is grid;
     
     In this function the list of neighboring atoms based on grid is generated.
     Input is grid;
     Output is atom.neighb;
  */



  for (Atom *atm = atom; atm != NULL; atm=atm->next) {
//    cout << "Molecule read_PDB " << atm->tag << " " << atm->resn << " " << atm->chain_id << " " << atm->resi 
//         << " " << atm->crd.x << " " << atm->crd.y << " " << atm->crd.z << " " << atm->het << " " << atm->model << endl;

    atm->delete_neighbor_list();
  }

  /* v gridu vsi atomi (tudi vsi modeli & hetero), ki so v mejah (0,0,0) (NUM_CELLS,NUM_CELLS,NUM_CELLS) - make_grid sam preveri!!! */
  for (Atom *atm = atom; atm != NULL; atm=atm->next) {
    grid->make_grid(atm);
  }

  /* samo za prvi model, query chain, in samo za ne-het atome naredimo listo neighborjev */
  for (Atom *atm = atom; atm != NULL; atm=atm->next) {
    if (atm->model == 1 && !atm->het && chain_id.find(atm->chain_id) != string::npos) grid->volume_slice_atom(atm, atm->r + distance + MAXR, distance); // modified OCT/22/2008 use 4*PROBE when testing! ( obsolete!! :)
  }

}

pair<Coor,Coor> Molecule::minmax_coor() {
  /*
    Dobimo najvecjo in najmanjso koordinato v query chain-u. Upostevamo
    samo prvi model (ne prileganih proteinov) in samo ne-hetero atome pravega (query)
    chain-a.
  */
  Coor minCrd (999999.0, 999999.0, 999999.0);
  Coor maxCrd (-999999.0, -999999.0, -999999.0);

  for (Atom *a = atom ; a!= NULL; a=a->next) {

#ifdef VERB
    cout << "x = " << a->crd.x << " model = " << a->model << " het = " << a->het << " chain = " << a->chain_id 
         << " chain_id = " << chain_id << endl;
#endif

    if (a->model == 1 && !a->het && chain_id.find(a->chain_id) != string::npos) {
      if (a->crd.x < minCrd.x) minCrd.x = a->crd.x;
      if (a->crd.y < minCrd.y) minCrd.y = a->crd.y;
      if (a->crd.z < minCrd.z) minCrd.z = a->crd.z;

      if (a->crd.x > maxCrd.x) maxCrd.x = a->crd.x;
      if (a->crd.y > maxCrd.y) maxCrd.y = a->crd.y;
      if (a->crd.z > maxCrd.z) maxCrd.z = a->crd.z;
      
    }
  }

  return make_pair(minCrd, maxCrd);
}

void Molecule::mark_motif_atoms(Motif *mf) {
  /*
    Oznacimo atome, ki ne pripadajo motivu v mf, z visited = true. Tako bo povrsina izracunana samo za
    atome motiva. Povrsina, to so probe centri.
   */
  for (Atom *a = atom; a != NULL; a=a->next) {
    if (!mf->is_in_motif(ChResi(a->chain_id, a->resi))) {
      a->visited = true;
    }
  }
}


void Molecule::all_triples(Probe *&probe) {
  Coor center_1, center_2;
  Atom *atm1, *atm2, *atm3;
  Probe *p;
  int p_size = 0;
  EElement *tmp3, *tmp4;
  for (atm1 = atom; atm1 != NULL; atm1=atm1->next) {
//    if (!atm1->neighb)
//      cout << "atom nima neighborjev" << endl;
//    else
//      cout << "neighborji so" << endl;
//    cout << "atm->num = " << atm1->num << endl;
    if (!atm1->visited) {
      atm1->visited = true;
      for (EElement *tmp = atm1->neighb; tmp != NULL; tmp = (EElement*) tmp->next) {
        atm2 = (Atom*) tmp->item;
        if (!atm2->visited) {
          for (EElement *tmp2 = (EElement*) tmp->next; tmp2 != NULL; tmp2 = (EElement*) tmp2->next) {
            atm3 = (Atom*) tmp2->item;
            if (!atm3->visited)
              if (atm3->r + 2*PROBE + atm2->r - dist(atm3->crd, atm2->crd) > EPS) {
                if (probe_center(atm1, atm2, atm3, center_1, center_2)) {
//                atm1->print_atom();
//                atm2->print_atom();
//                atm3->print_atom();
//                cout << " center_1 " << center_1.x << " " << center_1.y << " " << center_1.z << endl;
//                cout << " center_2 " << center_2.x << " " << center_2.y << " " << center_2.z << endl;
                
                  for (tmp3 = atm2->neighb; tmp3 != NULL; tmp3 = (EElement*) tmp3->next) {
                    if (tmp3->item == atm3) {
                      break;
                    }
                  }
                  for (tmp4 = atm3->neighb; tmp4 != NULL; tmp4 = (EElement*) tmp4->next) {
                    if (tmp4->item == atm2) {
                      break;
                    }
                  }
                  if (!atm1->overlap(center_1)) {
                    p = new Probe(p_size++);
                    p->crd = center_1;
                    p->next = probe;
#ifdef CILE
                    p->r = PROBE; // dolocimo probe radij
#endif
                    probe = p;

//                  cout << "tmp->size = " << tmp->size << endl;
//                  cout << "tmp2->size = " << tmp2->size << endl;
//                  cout << "tmp3->size = " << tmp2->size << endl;
//                  cout << "tmp4->size = " << tmp2->size << endl;

                    tmp->probe[tmp->size] = p;
                    tmp->atm[tmp->size++] = atm3;
                    
                    if (atm2->num < atm3->num) {
                      tmp3->probe[tmp3->size] = p;
                      tmp3->atm[tmp3->size++] = atm1;
                    }
                    else {
                      tmp4->probe[tmp4->size] = p;
                      tmp4->atm[tmp4->size++] = atm1;
                    }
                    
                    tmp2->probe[tmp2->size] = p;
                    tmp2->atm[tmp2->size++] = atm2;
                  }
                  if (!atm1->overlap(center_2)) {
                    p = new Probe(p_size++);
                    p->crd = center_2;
                    p->next = probe;
                    probe = p;
                    
//                  cout << "2 tmp->size = " << tmp->size << endl;
//                  cout << "2 tmp2->size = " << tmp2->size << endl;
//                  cout << "2 tmp3->size = " << tmp2->size << endl;
//                  cout << "2 tmp4->size = " << tmp2->size << endl;

                    tmp->probe[tmp->size] = p;
                    tmp->atm[tmp->size++] = atm3;
                  
                    if (atm2->num < atm3->num) {
                      tmp3->probe[tmp3->size] = p;
                      tmp3->atm[tmp3->size++] = atm1;
                    }
                    else {
                      tmp4->probe[tmp4->size] = p;
                      tmp4->atm[tmp4->size++] = atm1;
                    }
                    tmp2->probe[tmp2->size] = p;
                    tmp2->atm[tmp2->size++] = atm2;
                  }
                }
              }
          }
        }
        moves_from_centers(atm1, atm2, tmp);
      }
    }
    
  }
}

void Molecule::moves_from_centers(Atom *atm1, Atom *atm2, EElement *edge) {
  /*
    We have an edge with probe centers around it. We have to know which pairs of probe centers are on accessible from 
    one anoher.
    Output: In the array probe->center[i].move, we write a possible move of each probe->center[i] (there can be one for each edge)
  */
  if (edge->size % 2 != 0) { 
#ifdef VERB
    cout << "Warning (MOLECULE) : edge->size is odd " << edge->size << endl; 
#endif
    edge->size--; 
  }
  if (edge->size > 0) {
    Coor at21 = atm2->crd - atm1->crd;
    Coor X = atm1->crd + project(edge->probe[0]->crd - atm1->crd, at21);
    Coor center_new;
    double fi[10], tmp_fi;
    Coor vec[10], tmp_vec;
    double r;
    Coor x_vec, y_vec;

#ifdef VERB
    cout << "edge->size = " << edge->size << endl;
    cout << "edge->probe[0]->crd = (" << edge->probe[0]->crd.x << "," << edge->probe[0]->crd.y << "," << edge->probe[0]->crd.z << ")" << endl;
#endif    

    vec[0] = edge->probe[0]->crd - X;
    fi[0] = 0;
    for (int i = 1; i < edge->size; i++) {
      vec[i] = edge->probe[i]->crd - X;
      fi[i] = angle(vec[0], vec[i], at21);

    }
    x_vec = vec[0];
    y_vec =  at21 % x_vec;
    r = sqrt(x_vec * x_vec);

    x_vec = norm(x_vec);
    y_vec = norm(y_vec);
#ifdef VERB
    cout << "x_vec.x = " << x_vec.x << "x_vec.y = " << x_vec.y << "x_vec.z = " << x_vec.z << endl;
    cout << "y_vec.x = " << y_vec.x << "y_vec.y = " << y_vec.y << "y_vec.z = " << y_vec.z << endl;
#endif

    center_new = rotate_vector(x_vec, y_vec, X, 0.001, r);
    
    if (dist(center_new, edge->atm[0]->crd) > dist(edge->probe[0]->crd, edge->atm[0]->crd)) {
      for (int i = 1; i < edge->size; i++)
        for (int j = i + 1; j < edge->size; j++)
          if (fi[j] < fi[i]) {
            edge->exchange(i, j);
            tmp_fi = fi[i]; fi[i] = fi[j]; fi[j] = tmp_fi;
            tmp_vec = vec[i]; vec[i] = vec[j]; vec[j] = tmp_vec;
          }
    }
    else {
      fi[0] = 4 * acos(0);
      for (int i = 1; i < edge->size; i++)
        for (int j = i + 1; j < edge->size; j++)
          if (fi[j] > fi[i]) {
            edge->exchange(i, j);
            tmp_fi = fi[i]; fi[i] = fi[j]; fi[j] = tmp_fi;
            tmp_vec = vec[i]; vec[i] = vec[j]; vec[j] = tmp_vec;
          }
    }    
    
    for (int i = 0; i < edge->size; i = i + 2) {
      edge->probe[i]->move[edge->probe[i]->size_move++] = edge->probe[i+1];
      if (edge->probe[i]->size_move > 3) { 
	throw Err("Error (MOLECULE) : probe->size_move > 3 (=" + to_string(edge->probe[i]->size_move) + ")\n", 4);
      }
      edge->probe[i+1]->move[edge->probe[i+1]->size_move++] = edge->probe[i];
      if (edge->probe[i]->size_move > 3) { 
	throw Err("Error (MOLECULE) : probe->size_move > 3 (=" + to_string(edge->probe[i]->size_move) + ")\n", 4);
      }
    }
  }
}

void Molecule::expand(Probe *probe, int color) {
  if (probe->color == 0) {
    probe->color = color;
    surf_size_color[color]++;
    for (int i = 0; i < probe->size_move; i++) 
        expand(probe->move[i], color);
  }
}

void Molecule::enumerate_clefts(Probe *probe) {
  /*
    Here we enumerate the cavities. We give the probe centers in probe->center distinct colors indicating to which cavity 
    they belong to. We sum the number of the probe centers of each color in the surf_size_color[color]; Then, we search 
    for the color with maximum number of probe atoms == the largest surface and print this color in the lcolor variable.
  */
  int color = 1, max;
  Probe *p = probe;
  while (p != NULL) {
    if (p->color == 0) {
      surf_size_color[color] = 0;
      expand(p, color);
      color++;
    }
    p = p->next;
  }
  max = 0;
  for (int i = 1; i < color; i++) {
    if (surf_size_color[i] > max) { lcolor = i; max = surf_size_color[i]; }
    //    cout << i << "-th color size = " << surf_size_color[i] << endl;
  }
}






#ifdef CILE
void Molecule::surface_atoms_cile(Grid *grid, Probe *probe) {
  /*
     V atom.colors zapisemo enko pri vsakem atomu, ki je v kontaktu s probe-om, ki pripada danemu patch-u 
     s pribliznim radijem CILER.
  */

  // inicializiramo atom->colors
  for (Atom *a = this->atom; a != NULL; a = a->next) a->colors.size = 0;

  // gremo po probe-ih, ki pripadajo patchu
  Probe *p = probe;
  while (p != NULL) {
    if (p->color && p->dist < CILER) {
      // zapisemo v atom->colors barvo patch-a (=1)
      grid->volume_slice_surf(p, PROBE + MAXR + SURF, PROBE + SURF);
    }
    p = p->next;
  }
}

void Molecule::dijkstra(set<Probe*> Graph, Probe *source) {
  /*
    Dijkstra algoritem za iskanje najkrajsih poti. Izracunamo najkrajso razdaljo od tocke source do vseh tock v grafu,
    rezultati so v (Probe*)->dist za vsako tocko.
  */
  for (set<Probe*>::iterator it = Graph.begin(); it != Graph.end(); it++) { // Initializations
    (*it)->dist = XXL; // Unknown distance function from source to v
    (*it)->previous = NULL; // Previous node in optimal path from source
    (*it)->inQ = true; // Oznacimo odstranjene vertekse
  }
  source->dist = 0.0; // Distance from source to source
  set<Probe*> Q;
  Q.insert(Graph.begin(), Graph.end());
  // All nodes in the graph are unoptimized - thus are in Q
  while (!Q.empty()) { // The main loop
    set<Probe*>::iterator uit = min_element(Q.begin(), Q.end(), by_probe_dist);
    Probe *u = *uit;
    if (u->dist == XXL) {
      break ;                        // all remaining vertices are inaccessible from source
    }
    Q.erase(uit);
    u->inQ = false;
    for (int i = 0; i < u->size_move; i++)  { // where v has not yet been removed from Q.
      Probe *v = u->move[i];
      if (v->inQ) { // preverimo, ce nismo ze odstranili v-ja iz Q-ja
        float alt = u->dist + dist(u->crd, v->crd);
        if (alt < v->dist) {  // Relax (u,v,a)
          v->dist = alt ;
          v->previous = u;
        }
      }
    }
  }
}

void Molecule::expand_cile(Probe *pc, Probe *pd) {
  /*
    Gremo rekurzivno iz zacetnega pc-ja po vseh bliznjih probe centrih. Tiste, ki so preblizu zacetnemu pc-ju oznacimo kot
    "too close".
  */
  if (pd->color == 0 && dist(pc->crd, pd->crd) < CILER) {

    // oznacimo too close
    if (dist(pc->crd, pd->crd) < 2.0) {
      pd->too_close = true;
    }
    pd->color = 1;
    for (int i = 0; i < pd->size_move; i++) 
        expand_cile(pc, pd->move[i]);
  }
}

void Molecule::patches_cile(Grid *grid, Probe *probe) {
  /*
    Izpisemo vse patche na povrsini (v obliki probe centrov), ki imajo manhattan radius < CILER od 
    sredinskega probe centra. Manhattan distance je misljen kot nakrajsa pot med sredinskim probe
    centrom in obravnavanim probe centrom, tako da potujemo samo po ravnih crtah med probe centri.
  */

  // parametri za nastavljat
  const int atom_no = 610;  // stevilka centralnega atoma za patch, ki ga hocemo
  const float d = 3.0;  // izberemo tisti probe, ki je manj kot d od centralnega atoma

  int patch_no = 0;

  // inicializiramo
  for (Probe *t = probe; t != NULL; t = t->next) {
    t->too_close = false;
  }

  // dodatek, da izberemo samo en patch, ki je centriran na atomu dolocene stevilke
  Atom *ac=NULL;
  for (Atom *a = this->atom; a != NULL; a = a->next) {
    if ( a->num == atom_no) ac = a;
  }
  //////////////////////////

  // glavna zanka
  Probe *p = probe;
  while (p != NULL) {
    
    // preskocimo, ce ta p ni blizu centralnega atoma
    if (dist(p->crd, ac->crd) > d) { p=p->next; continue; }
    /////////////////////////////

    // da ne delamo patchev preblizu skupaj
    if (!p->too_close) {

//      // dolocimo neighb (atome) za p
//      grid->volume_slice_probe(p, PROBE + MAXR + SURF, PROBE + SURF);
//      
//      // ali smo katerega od sosednjih atomov ze obiskali?
//      cout << "--------------------- bliznji atomi od p-ja stevilka " << p->num << " -----------------" << endl;
//      for (EElement *pn = p->neighb; pn != NULL; pn = (EElement*) pn->next) {
//        ((Atom*) pn->item)->pdb();
//      }
//      
//      // ce ja, potem zaenkrat ne naredimo nic (vzamemo toliko patchev kolikor je probe-ov)
      
      // gremo po vseh moznih povezavah med p-ji, tako da ostanemo na isti SAS
      if (p->color == 0) {
        expand_cile(p, p);
      }
      
      // zapisemo najdeni graf v obliko, ki jo zna prebrati Dijkstra algoritem
      set<Probe*> Graph;
      for (Probe *t = probe; t != NULL; t = t->next) {
        if (t->color) {
          Graph.insert(t);
        }
      }
      
      // pozenemo dijkstra algoritem
      dijkstra(Graph, p);
      
      // poiscemo atome, ki so bliznji probe-om v patchu, in jim v color zapisemo 1
      surface_atoms_cile(grid, probe);
      
      // izpisemo patch kot PDB zapis probe-ov
      cout << "PATCH " << patch_no << " PROBES ";
      for (Probe *t = probe; t != NULL; t = t->next) {
        if (t->color) {
          cout << t->num << " ";
//          t->pdb();
        }
      }
      
      // izpisemo patch kot PDB zapis atomov
      cout << endl << "PATCH " << patch_no << " ATOMS ";
      for (Atom *a = this->atom; a != NULL; a = a->next) {
        if (a->colors.size > 0) {
          cout << a->num << " ";
//          a->pdb();
        }
      }
      // izpisemo patch kot PDB zapis resi-jev, lahko se podvajajo
      cout << endl << "PATCH " << patch_no << " RESI ";
      for (Atom *a = this->atom; a != NULL; a = a->next) {
        if (a->colors.size > 0) {
          cout << a->resi << " ";
//          a->pdb();
        }
      }

      cout << endl;
      
      patch_no++;

      // vsakic znova resetiramo color, ker moramo expand-at na vsakem probe centru
      for (Probe *t = probe; t != NULL; t = t->next) t->color = 0;
    }
    p = p->next;
  }
}
#endif // CILE







//void Molecule::output(bool _motif) {
void Molecule::output() {
  cout << "MOL> " << lcolor << endl;
  for (Atom *atm = atom; atm != NULL; atm = atm->next) {
    atm->print_atom();
  }
}

//void Molecule::output() {
//  cout << "MOL> " << lcolor << endl;
//  for (Atom *atm = atom; atm != NULL; atm = atm->next) {
//    if (atm->colors.size > 0)
//      atm->print_atom();
//  }
//}



void Molecule::read_plain_pdb(bool _firstModel) {
  /*
    Preberemo PDB file kot tekstovni file. Ce je _firstModel = true, potem preberemo 
    le prvi model (do prvega END), v nasprotnem primeru, beremo do konca (EOF).
   
   */

  ifstream pdb ((add_end_slash(INDIR) + name).c_str(), ios::binary);
  const size_t line_width = 80;

  if (! pdb.is_open()) { 
    throw Err("Error (MOLECULE) : Cannot open PDB file for plain reading " + name + ".", 5); 
  }


  // dobimo stevilo vrstic v file-u
  pdb.seekg(0, ios::end);
  size_t length = pdb.tellg();
  pdb.seekg(0, ios::beg);
  
#ifdef VERB
  cout << "Molecule:read_plain_pdb reserved space = " << (size_t) (length / line_width * 1.1) << endl;
#endif
  plain_pdb.reserve((size_t) (length / line_width * 1.1));

#ifdef VERB
  size_t c = 0;
#endif

  while (! pdb.eof() ) {
    string s;
    getline (pdb, s);

#ifdef VERB
    cout << "Molecule:read_plain_pdb capacity of s = " << s.capacity() << endl;
#endif
    /* appendaj toliko blank znakov, da bo dolzina stringa line_width (pride v postev pri PDB-fileih, ki nimajo beta faktorjev */
    if (line_width > s.size()) {
      s.append(line_width - s.size(), ' ');
    }

    if (!pdb.eof()) {
      plain_pdb.push_back(s);
    }

    if (_firstModel && plain_pdb.back().compare(0, 3, "END") == 0) break; // ne beremo vec po END

#ifdef VERB
    if (c != plain_pdb.capacity()) {
      c = plain_pdb.capacity();
      cout << "Molecule:read_plain_pdb capacity = " << plain_pdb.capacity() << endl;
    }
#endif

  }

  pdb.close();

}


void Molecule::output_iface_pdb() {
  /*
    The complex protein pdb file is marked, so that the residues that have been found  
    to belong to the actual interface have beta factors set to 1.00.
  */

  string name = add_end_slash(OUTDIR) + pdb_id + ".iface.pdb";

  ofstream f(name.c_str());

  for (unsigned int j = 0; j < plain_pdb.size(); j++) {
    if (plain_pdb[j].compare(0, 6, "HETATM") == 0) {

      plain_pdb[j].replace(60, 6, "  0.00");
    }
    else if (plain_pdb[j].compare(0, 4, "ATOM") == 0) {
      
      plain_pdb[j].replace(60, 6, "  0.00");
      
      int resi = atoi( plain_pdb[j].substr(22,4).c_str() );
      char c_id = plain_pdb[j].at(21);
 
      if ( imap.find(make_pair(c_id, resi)) != imap.end() )
        plain_pdb[j].replace(60, 6, "  1.00" );
      
    }
    f << plain_pdb[j] << endl;
    
  }
  f.close();
  
}

void Molecule::gen_cons_biounit() {
  /* 
     Input je navaden PDB file, prebrali smo samo prvi model in query chain_id-je. 
     Output so bioloski uniti (vsak v svojem file-u), tako da rotiramo in transliramo query chain id-je, ki so v dolocenem biounitu, 
     vsako rotacijsko matriko (BIOMT) zapisemo v svoj MODEL. Naredimo vse BIOMOLECULE, ki vsebujejo vsaj
     enega od obravnavanih chainov na query proteinu. 
     Za NMR strukture ne zgeneriramo; ce pozenemo za output_cons_pdb, potem imamo tudi stopnje ohranjenosti v beta faktorjih.
  */
  
  Kabsch *k = new Kabsch();
  
  
  /* gremo po vseh bioloskih unitih */
  for (map<pair<int,int>, BioUnit*>::iterator it = biounit.begin(); it != biounit.end(); it++) {
    
//    int model = it->first.first;
    int bio_num = it->first.second;
    BioUnit *bu = it->second;
    int new_model = 1;
    
    // naredimo biounit samo, ce se ta biounit nanasa na vsaj en chain_id iz this.chain_id 
    /* query this.chain_id je AB in bu->chain_id = {A,B,C} */
    if (bu->najdi(chain_id)) {
      // zapisemo ga v file
//      string name = add_end_slash(OUTDIR) + pdb_id + ".bu" + to_string(bio_num) + ".cons.pdb";
      string name = add_end_slash(OUTDIR) + pdb_id + chain_id + ".bu" + to_string(bio_num) + ".cons.pdb";

      bu->filename = name;

      ofstream f(name.c_str());
      
      /* gremo po vseh matrikah in rotiramo koordinate - za vsako matriko naredimo en MODEL - ENDMDL odsek */
      for (map<int, gsl_matrix*>::iterator it2 = bu->U.begin(); it2 != bu->U.end(); it2++) {
        
        int rota_num = it2->first;
        
        char mdl[SMALL];
        sprintf(mdl,"%-6s%8d%10s", "MODEL", new_model, (pdb_id + chain_id).c_str() );
        new_model++;
        
        bu->rota_pdb.push_back(mdl);
#ifdef VERB
        cout << bu->rota_pdb.back() << endl;
#endif      
        
        gsl_matrix *U = it2->second;
        gsl_vector *t = bu->t[rota_num];
        
#ifdef VERB
        cout << "generate_biounit() MATRIX " << it2->first << endl;
        cout << "generate_biounit() MATRIX " << gsl_matrix_get(U, 0, 0) << " " 
             << gsl_matrix_get(U, 0, 1) << " " << gsl_matrix_get(U, 0, 2) << " " << endl;
        cout << "generate_biounit() MATRIX " << gsl_matrix_get(U, 1, 0) << " " 
             << gsl_matrix_get(U, 1, 1) << " " << gsl_matrix_get(U, 1, 2) << " " << endl;
        cout << "generate_biounit() MATRIX " << gsl_matrix_get(U, 2, 0) << " " 
             << gsl_matrix_get(U, 2, 1) << " " << gsl_matrix_get(U, 2, 2) << " " << endl;
        cout << "generate_biounit() VECTOR " << gsl_vector_get(t, 0) << " " 
             << gsl_vector_get(t, 1) << " " << gsl_vector_get(t, 2) << " " << endl;
#endif      
        

        /* gremo po vseh atomih proteina, jih rotiramo in zapisemo v biounit->rota_pdb; plain_pdb smo prebrali z modifierjem _firstModel */
        for (unsigned int i = 0; i < plain_pdb.size(); i++) {
          
          /* ce se biomolecule chain in trenutno brani chain ujemata, oziroma, ce v bu->chain sploh ni nobenih
             chainov (NMR struktura), potem rotiramo koordinate */
          if (plain_pdb[i].compare(0, 4, "ATOM") == 0 || plain_pdb[i].compare(0, 6, "HETATM") == 0) {
            if (bu->is_in_rota_num(rota_num, plain_pdb[i].at(21))) {
//            if (bu->chain_id.find(plain_pdb[i].at(21)) != bu->chain_id.end()) {
              
              Coor c;
              c.x = atof( plain_pdb[i].substr(30,8).c_str() );
              c.y = atof( plain_pdb[i].substr(38,8).c_str() );
              c.z = atof( plain_pdb[i].substr(46,8).c_str() );
              
              /* rotiramo in transliramo */
              Coor crot = k->rotate_vector(c, U, t);
              
              char temp[SMALL];
              
              sprintf(temp,"%8.3f%8.3f%8.3f", crot.x, crot.y, crot.z);
              
              /* v rota_pdb kopiramo rotirane koordinate in vse ostalo */
              bu->rota_pdb.push_back( plain_pdb[i] );
              bu->rota_pdb.back().replace(30, 24, temp);
              
#ifdef VERB
              cout << bu->rota_pdb.back() << endl;
#endif      
              
            }
          }
	  /* zapisemo vse podatke, da lahko 1) identificiramo MODRES 2) SITE in njihov opis v REMARK 800 3) vse razlicne BIOMOLECULE */
	  else if (plain_pdb[i].compare(0, 6, "MODRES") == 0 || 
		   plain_pdb[i].compare(0, 6, "SITE  ") == 0 ||
		   plain_pdb[i].compare(0, 5, "LINK ") == 0 ||
//		   plain_pdb[i].compare(0, 3, "TER") == 0 ||
		   plain_pdb[i].compare(0, 10, "REMARK 800") == 0 ||
		   plain_pdb[i].compare(0, 10, "REMARK 350") == 0 ||
		   plain_pdb[i].compare(0, 6, "HETNAM") == 0) {
	    
            bu->rota_pdb.push_back( plain_pdb[i] );
//	    f << plain_pdb[i] << endl;
	  }
        }
        
        bu->rota_pdb.push_back("ENDMDL");


#ifdef VERB
        cout << bu->rota_pdb.back() << endl;
#endif      

      }
      bu->rota_pdb.push_back("END");
      // zapisemo file z biological unitom ( v beta faktorjih je ohranjenost, ce smo prej pognali output_cons_pdb )
      for (vector<string>::iterator it3 = bu->rota_pdb.begin(); it3 != bu->rota_pdb.end(); it3++)
        f << *it3 << endl;
      
      f.close();


    }


  }
  
  
  delete k;
  
}


void Molecule::remark_pdb(Item *one_clq) {
  /*
    Spremeni REMARK zapis v plain_pdb-ju tako, da bosta, ceprav je protein rotiran, 
    nova rotacijska matrika in translacijski vektor dala bioloski unit.
    Glej tudi zapiske, dne 10. marca 2010 :-)
    
    Situacija: Ta funkcija se uporablja z modifikatorjem --super, 
    po kateri je drugi protein 
    rotiran in transliran tako, da je superimposan na query (prvi protein). Matrike
    za zgenerirat bioloske unite drugega proteina so zato neveljavne in jih je treba 
    na novo izracunat, in sicer tako, da ko pomnozis koordinate atomov proteina s temi 
    novimi matrikami, dobis bioloski unit superimposanega proteina.

  */
  
  const int model = 1;
  int bio_num = 0;
  gsl_matrix *UR = gsl_matrix_alloc(3,3);
  gsl_vector *TR = gsl_vector_alloc(3);

  Output *o = new Output();

  for (size_t i = 0; i < plain_pdb.size(); i++) {


    if (plain_pdb[i].compare(0, 23, "REMARK 350 BIOMOLECULE:") == 0) {
      bio_num++;
    }

    /* zapisemo v plain_pdb nove rotacijske matrike U in translacijske vektorje t za ta biounit*/
    /* REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 */
    if (plain_pdb[i].compare(0, 18, "REMARK 350   BIOMT") == 0) {


      int biomt_num = atoi(plain_pdb[i].substr(18,1).c_str());

      int rota_num =  atoi(plain_pdb[i].substr(19,4).c_str());
        
      pair <int, int> bio_id;

      bio_id = make_pair(model, bio_num);
      
      /* izracunamo novo rotacijsko matriko in vektor */
      if (biomt_num == 1) {


        gsl_matrix *U1 = biounit[bio_id]->U[rota_num];
        gsl_vector *T1 = biounit[bio_id]->t[rota_num];


        gsl_matrix *U2 = one_clq->U;
        gsl_vector *T2 = one_clq->t;

#ifdef VERB
        cout << "U1" << endl;
        o->out_rota(U1);
        cout << "T1" << endl;
        o->out_trans(T1);

        cout << "U2" << endl;
        o->out_rota(U2);
        cout << "T2" << endl;
        o->out_trans(T2);
#endif

        /* nova matrika UR = U2(transposed) * U1 * U2 */
        gsl_matrix *URR = gsl_matrix_alloc(3,3);
        gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                        1.0, U2, U1,
                        0.0, URR);

        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, URR, U2,
                        0.0, UR);

        gsl_matrix_free(URR);


        /* nov vektor TR = U2(transposed) * (U1 * T2 + T1 - T2) */
        gsl_vector *TRR = gsl_vector_alloc(3);
        gsl_blas_dgemv(CblasNoTrans, 1, U1, T2, 0, TRR);

        gsl_vector_add(TRR, T1);

        gsl_vector_sub(TRR, T2);

        gsl_blas_dgemv(CblasTrans, 1, U2, TRR, 0, TR);

        gsl_vector_free(TRR);


#ifdef VERB
        cout << "UR" << endl;
        o->out_rota(UR);

        cout << "TR" << endl;
        o->out_trans(TR);
#endif
      }

      /* zapisemo nove rotacijske matrike in vektorje v plain_pdb */
      double URa = gsl_matrix_get(UR, biomt_num - 1, 0);
      double URb = gsl_matrix_get(UR, biomt_num - 1, 1);
      double URc = gsl_matrix_get(UR, biomt_num - 1, 2);
      
      double TRx = gsl_vector_get(TR, biomt_num - 1);
      
      stringstream s; 
      s.str("");s << setw(10) << setfill(' ') << fixed << setprecision(6) << URa;
      plain_pdb[i].replace(23, 10, s.str() );
      s.str("");s << setw(10) << setfill(' ') << fixed << setprecision(6) << URb;
      plain_pdb[i].replace(33, 10, s.str() );
      s.str("");s << setw(10) << setfill(' ') << fixed << setprecision(6) << URc;
      plain_pdb[i].replace(43, 10, s.str() );
      s.str("");s << setw(15) << setfill(' ') << fixed << setprecision(5) << TRx;
      plain_pdb[i].replace(53, 15, s.str() );
    }
      
  }

  gsl_matrix_free(UR);
  gsl_vector_free(TR);


  delete o;
}

string Molecule::remark_scores(Item *one_clq) {
  /* 
     izpisemo score kot REMARK
   */
  string s;
  s += "REMARK PROBIS ALIGNED_VERTICES " + to_string(one_clq->vert.size()) + "\n";
  s += "REMARK PROBIS E_VALUE " + to_string(one_clq->blosum_score) + "\n";
  s += "REMARK PROBIS RMSD " + to_string(one_clq->calpha_rmsd) + "\n";
  s += "REMARK PROBIS SVA " + to_string(one_clq->surf_vector_angle) + "\n";
  s += "REMARK PROBIS Z_SCORE " + to_string(z_score(one_clq->cluster_score)) + "\n";
  s += "REMARK PROBIS ALIGNMENT_SCORE " + to_string(one_clq->cluster_score) + "\n";
  return s;
}

void Molecule::output_rotated_asym(Molecule *m, Item *one_clq, int alignment_no) {
  /*
    Rotiramo biounit alignanega proteina na predstojnika in izpisemo samo alignan protein, z nekaterimi
    zapisi iz header sekcije.
  */
  string name = add_end_slash(OUTDIR) + m->pdb_id + m->chain_id + "_" + pdb_id + chain_id + "." + to_string(alignment_no) + ".rota.pdb";
  ofstream f(name.c_str());
  Kabsch *k = new Kabsch();

  f << remark_scores(one_clq); // izpisemo score kot remakr
  
  for (unsigned int i = 0; i < plain_pdb.size(); i++) {
    if (plain_pdb[i].compare(0, 4, "ATOM") == 0 || plain_pdb[i].compare(0, 6, "HETATM") == 0) {
      
      Coor c;
      
      c.x = atof( plain_pdb[i].substr(30,8).c_str() );
      c.y = atof( plain_pdb[i].substr(38,8).c_str() );
      c.z = atof( plain_pdb[i].substr(46,8).c_str() );
      
      /* we have to use inverse rotation because we are rotating PROTEIN2 to PROTEIN1 */
      Coor crot = k->rotate_vector_inv(c, one_clq->U, one_clq->t);
      char temp[SMALL];
      sprintf(temp,"%8.3f%8.3f%8.3f", crot.x, crot.y, crot.z);
      string s = plain_pdb[i];

      s.replace(30, 24, temp);
      f << s << endl;
    }
    /* zapisemo vse podatke, da lahko 1) identificiramo MODRES 2) SITE in njihov opis v REMARK 800 3) vse razlicne BIOMOLECULE */
    else if (plain_pdb[i].compare(0, 6, "MODRES") == 0 || 
//             plain_pdb[i].compare(0, 5, "MODEL") == 0 ||
             plain_pdb[i].compare(0, 3, "END") == 0 ||
//             plain_pdb[i].compare(0, 6, "ENDMDL") == 0 ||
             plain_pdb[i].compare(0, 6, "SITE  ") == 0 ||
//             plain_pdb[i].compare(0, 3, "TER") == 0 ||
             plain_pdb[i].compare(0, 10, "REMARK 800") == 0 ||
             plain_pdb[i].compare(0, 10, "REMARK 350") == 0 ||
             plain_pdb[i].compare(0, 6, "HETNAM") == 0) {

      f << plain_pdb[i] << endl;
    }
  }
  
  f.close();

  delete k;

}


void Molecule::output_rotate_pdb(Molecule *m, Item *one_clq) {
  /*
    Rotate a protein and output it together with protein m.
  */
  string name = add_end_slash(OUTDIR) + m->pdb_id + m->chain_id + "_" + pdb_id + chain_id + "." + to_string(ALIGNMENT_NO) + ".rota.pdb";
  ofstream f(name.c_str());
  Kabsch *k = new Kabsch();
  char mdl[SMALL];
  
  sprintf(mdl,"%-6s%8d%10s", "MODEL", 1, (pdb_id + chain_id).c_str() );
  f << mdl <<endl;

  for (unsigned int i = 0; i < plain_pdb.size(); i++) {
    if (plain_pdb[i].compare(0, 4, "ATOM") == 0 || plain_pdb[i].compare(0, 6, "HETATM") == 0) {
      
      Coor c;
      
      c.x = atof( plain_pdb[i].substr(30,8).c_str() );
      c.y = atof( plain_pdb[i].substr(38,8).c_str() );
      c.z = atof( plain_pdb[i].substr(46,8).c_str() );
      
      /* we have to use inverse rotation because we are rotating PROTEIN2 to PROTEIN1 */
      Coor crot = k->rotate_vector_inv(c, one_clq->U, one_clq->t);
      
      char temp[SMALL];
      sprintf(temp,"%8.3f%8.3f%8.3f", crot.x, crot.y, crot.z);
      plain_pdb[i].replace(30, 24, temp);


      /* popravimo beta faktor pri residujeu, ki je del patcha na 1.00 */
      char c_id = plain_pdb[i].at(21);
      int resi = atoi( plain_pdb[i].substr(22, 4).c_str() );

      plain_pdb[i].replace(60, 6, "  0.00");

      /* oznacimo z bfactor 1.00 tiste aminokisline, ki so ohranjene (samo na pravem chain-u!) */
      if (find(one_clq->cons.begin(), one_clq->cons.end(), make_pair(ChResi(c_id, resi), 2) ) != one_clq->cons.end() )
        plain_pdb[i].replace(60, 6, "  1.00");

      f << plain_pdb[i] << endl;
    }
    /* zapisemo vse podatke, da lahko 1) identificiramo MODRES 2) SITE in njihov opis v REMARK 800 3) vse razlicne BIOMOLECULE */
    else if (plain_pdb[i].compare(0, 6, "MODRES") == 0 || 
             plain_pdb[i].compare(0, 6, "SITE  ") == 0 ||
             plain_pdb[i].compare(0, 6, "HELIX ") == 0 ||
             plain_pdb[i].compare(0, 6, "SHEET ") == 0 ||
             plain_pdb[i].compare(0, 5, "LINK ") == 0 ||
//             plain_pdb[i].compare(0, 10, "REMARK 800") == 0 ||
             plain_pdb[i].compare(0, 10, "REMARK 350") == 0) {
//             plain_pdb[i].compare(0, 10, "REMARK 350") == 0 ||
//             plain_pdb[i].compare(0, 6, "HETNAM") == 0) {

      f << plain_pdb[i] << endl;

    }

    
  }
  
  f << "TER" << endl;
  f << "ENDMDL" << endl;
  
  
  
  sprintf(mdl,"%-6s%8d%10s", "MODEL", 2, (m->pdb_id + m->chain_id).c_str() );
  f << mdl <<endl;
  
  /* output the fixed protein as is */
  for (unsigned int i = 0; i < m->plain_pdb.size(); i++) {

    if (m->plain_pdb[i].compare(0, 4, "ATOM") == 0 || m->plain_pdb[i].compare(0, 6, "HETATM") == 0) {

      /* popravimo beta faktor pri residueju, ki je del patcha na 1.00 */
      char c_id = m->plain_pdb[i].at(21);
      int resi = atoi( m->plain_pdb[i].substr(22, 4).c_str() );


      m->plain_pdb[i].replace(60, 6, "  0.00");

      
      /* oznacimo z bfactor 1.00 tiste aminokisline, ki so ohranjene (samo na pravem chain-u!) */
      if (find(one_clq->cons.begin(), one_clq->cons.end(), make_pair(ChResi(c_id, resi), 1) ) != one_clq->cons.end())
        m->plain_pdb[i].replace(60, 6, "  1.00");



      f << m->plain_pdb[i] << endl;
    }
    /* zapisemo vse podatke, da lahko 1) identificiramo MODRES 2) SITE in njihov opis v REMARK 800 3) vse razlicne BIOMOLECULE */
    else if (m->plain_pdb[i].compare(0, 6, "MODRES") == 0 || 
             m->plain_pdb[i].compare(0, 6, "SITE  ") == 0 ||
             m->plain_pdb[i].compare(0, 6, "HELIX ") == 0 ||
             m->plain_pdb[i].compare(0, 6, "SHEET ") == 0 ||
//             m->plain_pdb[i].compare(0, 10, "REMARK 800") == 0 ||
             m->plain_pdb[i].compare(0, 10, "REMARK 350") == 0) {
//             m->plain_pdb[i].compare(0, 10, "REMARK 350") == 0 ||
//             m->plain_pdb[i].compare(0, 10, "HETNAM") == 0) {
//             m->plain_pdb[i].compare(0, 3, "TER") == 0) {

      f << m->plain_pdb[i] << endl;

    }

  }


  f << "TER" << endl;
  f << "ENDMDL" << endl;
  f << "END" << endl;

  f.close();

  delete k;

}


void Molecule::output_cons_pdb() {
  /*
    The query protein pdb or pdbqt file is marked, so that the residues that have been found conserved 
    have beta factors set to their conservation level.
  */

  string name = add_end_slash(OUTDIR) + pdb_id + ".cons.pdb";
  
  ofstream f(name.c_str());

  for (unsigned int j = 0; j < plain_pdb.size(); j++) {
    if (plain_pdb[j].compare(0, 6, "HETATM") == 0) {

      plain_pdb[j].replace(60, 6, "  0.00");
    }
    else if (plain_pdb[j].compare(0, 4, "ATOM") == 0) {
      
      plain_pdb[j].replace(60, 6, "  0.00");
      
      char c_id = plain_pdb[j].at(21);
      int resi = atoi( plain_pdb[j].substr(22,4).c_str() );

      if (chain_id.find(c_id) != string::npos) {
        Residue r(c_id, resi);
        if ( sequence.find(&r) != sequence.end() )

          plain_pdb[j].replace(64, 1, to_string( (*sequence.find(&r))->cons ) );
      }        
    
    }

    
    f << plain_pdb[j] << endl;

  }
  f.close();
  
}



void Molecule::out_query() {
  /*
    Izpisemo query.json, v katerem je zaporedje query proteina, ohranjenost, in podatek, ali je fingeprint.
  */
  
  ofstream f((add_end_slash(OUTDIR) + "query.json").c_str());
  
  f << "[";
  for (set<Residue*, rescomp>::iterator it = sequence.begin(); it != sequence.end(); it++) {
    if (it != sequence.begin()) f << ",";
    f << "{";
    f << "\"resi\":" << (*it)->resi<< ",";
    f << "\"resn\":\"" << (*it)->resn<< "\",";
    f << "\"chain_id\":\"" << (*it)->chain_id<< "\",";
    f << "\"cons\":" << (*it)->cons<< ",";
    f << "\"fp\":" << (*it)->sig_id;
    f << "}";
  }
  f << "]";
  f.close();
}



void Molecule::interface(float interchain_distance) {

  imap.clear();

  Atom *atm2;

#ifdef VERB
  char line[100];
#endif
  for (Atom *atm = atom; atm != NULL; atm = atm->next) {
    for (EElement *tmp = atm->neighb; tmp != NULL; tmp = (EElement*) tmp->next) {
      atm2 = (Atom*) tmp->item;
      if (atm->chain_id != atm2->chain_id && 
          atm->r + atm2->r + interchain_distance > dist(atm->crd, atm2->crd)) {

//        imap[make_pair(atm->chain, atm->resi)] = atm2->chain;

        /* '@' oznacuje protein chain in je rezervirana (query ne sme biti '@') */
        char veriga = atm2->chain_id;
        if (isalpha(veriga)) veriga = '@';

        /* ker vsak residue lahko pripada vecim interfaceom,e.g, protein-protein in protein-small ligand hkrati */
        pair<multimap<pair<char,int>,char>::iterator, multimap<pair<char,int>,char>::iterator> ret = imap.equal_range(make_pair(atm->chain_id,atm->resi));

        /* vendar moramo paziti, saj ta imap ima sedaj zapis za vsak atom, ki si je blizu - prej le za vsak residue */
        bool dodaj = true;

        /* moramo po vseh razlicnih b.s.-jih katerim en residue lahko pripada */
        for (multimap<pair<char,int>,char>::iterator it = ret.first; it != ret.second; it++) {
          if (it->second == veriga) { dodaj = false; break; }
        }


        /* dodamo samo, ce imap zapis ze ne obstaja, npr. ce je (A, 123) -> B ze not, ga ne dodas se enkrat  */
        if (dodaj) imap.insert ( make_pair(pair<char,int>(atm->chain_id,atm->resi), veriga) );




#ifdef VERB
	char line[1000];
        sprintf(line, "INTE> %c%c%4c%4d", atm->chain_id, veriga, one_letter_code(atm->resn), atm->resi);
        cout << line << endl;
#endif
      }
    }
  }
}


void Molecule::interface_models(float interchain_distance) {
  /*
    Prebrali smo .lig.pdb file, v katerem je predstavniki in vsi prilegani podobnezi z ligandi vred.
    Sedaj pa izracunamo, kateri ligandi so v prostoru blizu povrsini proteina-predstavnika. Posebna
    pozornost gre SITE zapisom v PDB, ti so obravnavani tukaj, a posebej od ostalih ligandov. 
    MODRES, ki so spremenjene aminokisline, sploh ne preberemo (izlocimo) ze v read_PDB.

    NOTE 18/APR/2010 : Tukaj je itak en sam chain_id
  */

  
  Atom *atm2;

  int iResi = -999999;
  char tmp_chain_id = ' ';

  for (Atom *atm = atom; atm != NULL; atm = atm->next) {
    /* binding site-e racunamo samo na query proteinu (pravi chain && model=1 && !het) */
    if (atm->model == 1 && !atm->het && chain_id.find(atm->chain_id) != string::npos) {
      /* upostevamo tudi SITE zapise na query chain-u (predstavniku) */
      /* samo enkrat na residue dodamo */
      if (iResi != atm->resi || tmp_chain_id != atm->chain_id) {
        tmp_chain_id = atm->chain_id;
        iResi = atm->resi;
        map<Ligand, string>::iterator sit = site.find( Ligand(atm->model, atm->chain_id, atm->resi) );
        if (sit != site.end() ) {
          string siteName = sit->second;
          /* vstavimo SITE na query chain-u: atm->resi (=-999) oznacuje SITE zapis */
          imap_models.insert ( make_pair( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi), 
                                          Ligand(atm->model, modelName[atm->model], atm->chain_id, siteName, -999)) );
        }
      }
      
      for (EElement *tmp = atm->neighb; tmp != NULL; tmp = (EElement*) tmp->next) {
        atm2 = (Atom*) tmp->item;
        /* pri SITE-ih upostevamo samo prilegane verige - tiste, katerim smo minimizirali rmsd glede na predstojnika ! */
        if (modelName[atm2->model].size() == 0 || atm2->chain_id == modelName[atm2->model].at(4) ) {
	  
          /* ce bliznji atom pripada SITE-u (zaenkrat vsakemu residue-ju lahko pripada le en site), 
             potem je kriterij za oddaljenost strozji (<3.0A) */
          map<Ligand, string>::iterator sit = site.find( Ligand(atm2->model, atm2->chain_id, atm2->resi) );
          if (!atm2->het && sit != site.end() && dist(atm->crd, atm2->crd) < SITE_RMSD) {
	    string siteName = sit->second;
#ifdef VERB
            char line[100];
            sprintf(line, "INTE> %6d%6s%6s%6d%6c --> %6d%6s%6s%6d%6c", atm->model, modelName[atm->model].c_str(), atm->resn, atm->resi, atm->chain_id, atm2->model, modelName[atm2->model].c_str(), siteName.c_str(), atm2->resi, atm2->chain_id);
            cout << line << " " << dist(atm->crd, atm2->crd) << endl;
#endif
            /* preden dodamo v imap_models preverimo, ce odnos med tem query residue in tem site-om ze obstaja 
               (vsak query residue lahko pripada vec binding siteom) */
            pair<multimap<Ligand,Ligand>::iterator, multimap<Ligand,Ligand>::iterator> ret = 
              imap_models.equal_range( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi) );
            
            /* paziti moramo, saj iscemo residueje, ki so si blizu skupaj (podatke pa imamo za atome) */
            bool dodaj = true;
            
            /* moramo po vseh razlicnih b.s.-jih katerim en residue lahko pripada */
            for (multimap<Ligand,Ligand>::iterator it = ret.first; it != ret.second; it++) {
              Ligand r2 = it->second;
              if ( atm2->model == r2.model && siteName.compare(r2.resn) == 0) { dodaj = false; break; }
            }
            /* dodamo samo, ce imap_models zapis ze ne obstaja (razlikujemo med hetero in proteinskimi atomi) */
            if (dodaj) {
              imap_models.insert ( make_pair( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi), 
                                              Ligand(atm2->model, modelName[atm2->model], atm2->chain_id, siteName, -999)) ); // -999 = SITE zapis !!
            }
          }
        }
#ifdef VERB        
        cout << "Molecule::interface_models()  : " << modelName[atm2->model] << " " << atm2->resn << " " << atm2->resi << " " << atm2->chain_id << endl;
#endif
        /* sedaj obravnavamo se prave ligande s 3D koordinatami */
        if (
            /* ce sta razlicni stevilki modela (e.g. 1,3), ce sta razlicna chaina (e.g. A,B), ali ce je sosed hetero atom */
            (atm->model != atm2->model || atm->chain_id != atm2->chain_id || atm2->het) &&
	    
            /* pri iskanju bliznjih proteinov-ligandov ne smemo upostevati verige, ki smo jo prilegali (samo tiste v model=2)! */
//            (atm2->het || !(atm2->model == 2 && atm2->chain_id == modelName[atm2->model].at(4))) &&
            (atm2->het || modelName[atm2->model].size() == 0 || !(atm2->model == 2 && atm2->chain_id == modelName[atm2->model].at(4))) &&

            /* vsi atomi imajo radije (hetero imajo r(hetero)=CATOM=1.5) */
            atm->r + atm2->r + interchain_distance > dist(atm->crd, atm2->crd)
            
           ) 
          {

#ifdef VERB        
          cout << "Molecule::interface_models()  : inside " << atm2->het << " " << modelName[atm2->model] << " " << atm2->resn << " " << atm2->resi << " " << atm2->chain_id << endl;
#endif

                
            /* preveriti moramo, ali ta odnos med query residue in ligandom ze obstaja (vsak query residue lahko pripada vec binding siteom) */
            pair<multimap<Ligand,Ligand>::iterator, multimap<Ligand,Ligand>::iterator> ret = 
              imap_models.equal_range( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi) );
            
            
            /* vendar moramo paziti, saj iscemo residueje, ki so si blizu skupaj (delamo pa na ravni atomov) */
            bool dodaj = true;
            int resi = (atm2->het == 0 ? 1 : atm2->resi);

            /* moramo po vseh razlicnih b.s.-jih katerim en residue lahko pripada */
            for (multimap<Ligand, Ligand>::iterator it = ret.first; it != ret.second; it++) {
              
              Ligand r2 = it->second;
                /* razlikujemo po stevilki modela, stevilki aminokisline (proteinski ligand= 1) in oznaki chaina (e.g., ne sme biti dveh z model=1 resn=301 ch=A) */
                if ( r2.model == atm2->model && r2.resi == resi && r2.chain_id == atm2->chain_id ) { dodaj = false; break; }
            }
            
#ifdef VERB        
            cout << "Molecule::interface_models()  : " << atm2->resn << " " << atm2->het << " " << dodaj << endl;
#endif
            
            /* dodamo samo, ce imap_models zapis ze ne obstaja (razlikujemo med hetero in proteinskimi atomi) */
            if (dodaj) {
              
              /* nucleic */
              if (!atm2->het && nucleic.find(atm2->resn) != string::npos) {
//                cout << "nucleic" << endl;

                /* vsaka DNK/RNK dobi oznake resn=+++ in resi=1 */
                imap_models.insert ( make_pair( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi), 
                                                Ligand(atm2->model, modelName[atm2->model], atm2->chain_id, "+++", 1)) );

              }
              /* protein */
              else if (!atm2->het && amino.find(atm2->resn) != string::npos) { 
//                cout << "protein" << endl;

                /* vsaka aminokislina (protein) dobi oznake resn=@@@ in resi=1 */
                imap_models.insert ( make_pair( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi), 
                                                Ligand(atm2->model, modelName[atm2->model], atm2->chain_id, "@@@", 1)) );

              }
              /* hetero */
              else  // if (ions.find(atm2->resn) != string::npos)
                imap_models.insert ( make_pair( Ligand(atm->model, modelName[atm->model], atm->chain_id, atm->resn, atm->resi), 
                                                Ligand(atm2->model, modelName[atm2->model], atm2->chain_id, atm2->resn, atm2->resi)) );
              
            


#ifdef VERB
              char line[100];
              sprintf(line, "INTE> %6d%6s%6s%6d%6c --> %6d%6s%6s%6d%6c", atm->model, modelName[atm->model].c_str(), atm->resn, atm->resi, atm->chain_id, 
                      atm2->model, modelName[atm2->model].c_str(), atm2->resn, atm2->resi, atm2->chain_id);
              cout << line << endl;
#endif
            }


          }

      }
    }

  }
}

void Molecule::output_interface_models() {
  /*
    Izpisemo v taksnile obliki :


       PREDSTOJNIK                                            LIGAND IZ PRILEGANEGA
       ----------------------------------------------------------------------------------------------------
            #MODEL     PDB        RESN     #RESI     CHAIN    #MODEL     PDB        RESN     #RESI      CHAIN   
       ----------------------------------------------------------------------------------------------------
       01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
       LIG>      1     1i52A       LEU       227         A         5     1injA       HOH      1033         A
  */

  for (multimap<Ligand, Ligand>::iterator it = imap_models.begin(); it != imap_models.end(); it++) {
    
    Ligand r1 = it->first;
    Ligand r2 = it->second;
    
    
    printf("LIG>%7d%10s%10s%10d%10c%10d%10s%10s%10d%10c\n", r1.model, 
           r1.pdb_id.c_str(), r1.resn.c_str(), r1.resi, r1.chain_id, 
           r2.model, r2.pdb_id.c_str(), r2.resn.c_str(), r2.resi, r2.chain_id);
    
  }
  
}

void Molecule::mark_binding_site(Molecule *m, string bsite) {
  /*
    Oznacimo atome, ki ne pripadajo binding site-u, tako da jih ne upoštevamo pri izračunu povrsine. 
    Binding site je podan z identifikatorjem liganda, na primer:
    *.*.A
    ATP.305.A
  */

  vector<string> s = reg_split(bsite, ".", "", "");
  if (s.size() != 3)
    throw Err("Error (MOLECULE) : Wrong bsite format : " + bsite + ".", 3);
#ifdef VERB
  cout << "Molecule::mark_binding_sites bsite=" << bsite << " resn="<< s[0] << " resi=" << s[1] << " chain_id=" << s[2] << endl;
#endif

  string resn = s[0];
  string resi = s[1];
  string chain_id = s[2];
  set<ChResi> sele;

  // najdemo vzorec residue-jev med razlicnimi binding site-i in ga zapisemo v sele
  for (multimap<Ligand, Ligand>::iterator it = m->imap_models.begin(); it != m->imap_models.end(); it++) {

    Ligand r1 = it->first;
    Ligand r2 = it->second;

#ifdef VERB
    cout << "Molecule::mark_binding_sites resi=" << r2.resi << " chain_id=" << r2.chain_id << endl;
    cout << "r2.resn = [" << r2.resn << "] resn = [" << resn << "]" << endl;
#endif
    if ((resn == "*" || r2.resn.find(resn)!=string::npos) && (resi == "*" || r2.resi == atoi(resi.c_str())) && (chain_id == "*" || r2.chain_id == chain_id.at(0))) {
      sele.insert(ChResi(r1.chain_id, r1.resi));
    }
  }

#ifdef VERB
  for (set<ChResi>::iterator it = sele.begin(); it != sele.end(); it++) {
    cout << "Molecule::mark_binding_sites sele resi=" << it->second << " chain_id=" << it->first << endl;
  }
#endif

  // oznacimo atome, ki ne pripadajo binding site-u
  for (Atom *a = atom; a != NULL; a=a->next) {
    if (sele.find(ChResi(a->chain_id, a->resi)) == sele.end()) {
      a->visited = true;
    }
  }
}

void Molecule::fn_important_residues() {
  /*
    V FNIMPORT zapisu v PDB fileu kompleksa lahko podamo residuje, ki so pomembni za funkcijo proteina, a ne pripadajo
    nobenemu interfaceu. Na primer iz literature zvemo, da je nek loop kriticno pomemben za 
    zvijanje (folding) tega proteina.
    0........10....16
    FNIMPORT  A     1   
    FNIMPORT  A    21   

   */
  for (unsigned int i = 0; i < plain_pdb.size(); i++) {
    if (plain_pdb[i].compare(0, 8, "FNIMPORT") == 0) {
      
      int resi = atoi( plain_pdb[i].substr(11,6).c_str() );

      char c_id = plain_pdb[i].at(10);

//      imap[make_pair(c_id, resi)] = '*';
      imap.insert ( make_pair(pair<char,int>(c_id,resi), '*') );

    }
  }

}


void Molecule::output_surface() {
  /*
    Output all surface atoms colored with color.
  */
  char line[300];
  char tmp[20];
  int i = 0;
  //  for (int i = 0; i < size; i++) {
  for (Atom *atm = atom; atm != NULL; atm = atm->next) {
    if (atm->colors.is_visited(lcolor)) {
      sprintf(line, "HETATM                                                                                               ");
      sprintf(tmp, "%d", ++i);
      sprintf(&line[11 - strlen(tmp)], "%s                                                                           ", tmp);
      sprintf(&line[13], "%s                                                                             ", atm->tag);
      sprintf(&line[17], "%s                                                                            ", atm->resn);
      sprintf(tmp, "%d", atm->resi);
      sprintf(&line[26 - strlen(tmp)], "%s                                                                           ", tmp);
      sprintf(tmp, "%2.3f", atm->crd.x);
      sprintf(&line[38 - strlen(tmp)], "%s                                                                           ", tmp);
      sprintf(tmp, "%2.3f", atm->crd.y);
      sprintf(&line[46 - strlen(tmp)], "%s                                                                           ", tmp);
      sprintf(tmp, "%2.3f", atm->crd.z);
      sprintf(&line[54 - strlen(tmp)], "%s  1.00                                                                     ", tmp);
      sprintf(tmp, "%2.2f", (double) lcolor + 8.0);
      sprintf(&line[66 - strlen(tmp)], "%s                                                                           ", tmp);
      line[66] = '\0';
      cout << line << endl;
    }
  }
  cout << "END" << endl;
}

void Molecule::surface_atoms(Grid *grid, Probe *probe) {
  /*
    Find all atoms of the protein, which are less than SURF Angstroms (see defs.h) below the surface. 
    The distance begins at the radius of the probe center and ends at the surface of the protein atom.
    Output is in atom.colors. Surface atoms are classified as to which closed surface they belong to 
    (atom.colors).
    This function works much in the same way as volume_slice().
    
  //    (NOT TRUE!) Output: besides above, in probe neighb array are surface atoms less than SURF Angstroms from 
  //    this probe center. Need for coloring the surface atoms after finding a clique by their geometric 
  //    complementarity.
  */
  double max_dist;
  Probe *p = probe;
  while (p != NULL) {
    max_dist = (PROBE + SURF + MAXR);
    grid->volume_slice_surf(p, max_dist, PROBE + SURF);
    p = p->next;
  }
}


void Molecule::set_descriptors(Descriptor *&desc, bool _nobb) {
  /*
    We set the schmitt descriptors to atoms in the molecule (atom[]). We consider each aminoacid separately.
    At the end, we set the desc.s.w[] arrays of neighboring schmitt descriptors to each descriptor.

    Preverimo tudi, ali upostevamo backbone atome, ali ne: _nobb je true ali false.

  */
  Atom *acid_start, *atm;
  Descriptor *d = NULL;
  bool done_pi, done_al;
  int d_size = 0;
  atm = atom;
  cout << "Setting descriptors ..." << endl;
  while (atm != NULL) {
    acid_start = atm->acid_start;
    done_pi = false;
    done_al = false;
    while(atm != NULL && atm->acid_start == acid_start) {
      //      if (atm->colors.size > 0) {
      if (atm->colors.is_visited(lcolor)) {
        // DONORS
        if (strcmp(atm->tag, "NE") == 0 || strcmp(atm->tag, "NH1") == 0 || strcmp(atm->tag, "NH2") == 0 || //ARG
            strcmp(atm->tag, "ND2") == 0 || //ASN
            (strcmp(atm->tag, "NE2") == 0 && strcmp(atm->resn, "GLN") == 0) || //GLN
            strcmp(atm->tag, "NZ") == 0 || //LYS
            strcmp(atm->tag, "NE1") == 0 || //TRP
            (!_nobb && strcmp(atm->tag, "N") == 0 ) ) //PEP
//            strcmp(atm->tag, "NE1") == 0) //TRP
          {
            d = new Descriptor(d_size++);
            d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
            d->s->mnsp = 0x08;
            d->psurf = 0;
            //            d->colors.copy(atm->colors);
            d->atom = atm;
            d->next = desc;
            desc = d;
          }      
        // ACCEPTORS
        if (strcmp(atm->tag, "OD1") == 0 || //ASN
            strcmp(atm->tag, "OD1") == 0 || strcmp(atm->tag, "OD2") == 0 || //ASP
            strcmp(atm->tag, "OE1") == 0 || //GLN
            strcmp(atm->tag, "OE1") == 0 || strcmp(atm->tag, "OE2") == 0 || //GLU
            (!_nobb && strcmp(atm->tag, "O") == 0 ) ) //PEP
//            strcmp(atm->tag, "OE1") == 0 || strcmp(atm->tag, "OE2") == 0) //GLU
          {
            d = new Descriptor(d_size++);
            d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
            d->s->mnsp = 0x04;
            d->psurf = 0;
            //            d->colors.copy(atm->colors);
            d->atom = atm;
            d->next = desc;
            desc = d;
          }
        // ACCEPTORS - DONORS
        if ((strcmp(atm->tag, "ND1") == 0 || strcmp(atm->tag, "NE2") == 0) && strcmp(atm->resn, "HIS") == 0) //HIS
          {
            d = new Descriptor(d_size++);
            d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
            d->s->mnsp = 0x0c;
            d->psurf = 0;
            //            d->colors.copy(atm->colors);
            d->atom = atm;
            d->next = desc;
            desc = d;
          }
        if (strcmp(atm->tag, "OG") == 0 && strcmp(atm->resn, "SER") == 0) //SER
          {
            d = new Descriptor(d_size++);
            d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
            d->s->mnsp = 0x0c;
            d->psurf = 0;
            //            d->colors.copy(atm->colors);
            d->atom = atm;
            d->next = desc;
            desc = d;
          }
        if (strcmp(atm->tag, "OG1") == 0 && strcmp(atm->resn, "THR") == 0) //THR
          {
            d = new Descriptor(d_size++);
            d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
            d->s->mnsp = 0x0c;
            d->psurf = 0;
            //            d->colors.copy(atm->colors);
            d->atom = atm;
            d->next = desc;
            desc = d;
          }
        if (strcmp(atm->tag, "OH") == 0 && strcmp(atm->resn, "TYR") == 0) //TYR
          {
            d = new Descriptor(d_size++);
            d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
            d->s->mnsp = 0x0c;
            d->psurf = 0;
            //            d->colors.copy(atm->colors);
            d->atom = atm;
            d->next = desc;
            desc = d;
          }
        // ALIPHATIC
        if (!done_al)
          if ((strcmp(atm->tag, "CB")  == 0 && strcmp(atm->resn, "ALA") == 0) || //ALA
              (strcmp(atm->tag, "CG")  == 0 && strcmp(atm->resn, "ARG") == 0) || //ARG
              (strcmp(atm->tag, "SG")  == 0 && strcmp(atm->resn, "CYS") == 0) || //CYS
              (strcmp(atm->tag, "CB")  == 0 && strcmp(atm->resn, "ILE") == 0) || //ILE
              (strcmp(atm->tag, "CB")  == 0 && strcmp(atm->resn, "LEU") == 0) || //LEU
              (strcmp(atm->tag, "SD")  == 0 && strcmp(atm->resn, "MET") == 0) || //MET
              (strcmp(atm->tag, "CG")  == 0 && strcmp(atm->resn, "PRO") == 0) || //PRO
              (strcmp(atm->tag, "CG2") == 0 && strcmp(atm->resn, "THR") == 0) || //THR
              (strcmp(atm->tag, "CB")  == 0 && strcmp(atm->resn, "VAL") == 0)) //VAL
            {                                     
              d = new Descriptor(d_size++);
              d->crd = atm->crd; //kopiranje vektorjev (operator= v geo.cc)
              d->s->mnsp = 0x01;
              d->psurf = 0;
              //              d->colors.copy(atm->colors);
              d->atom = atm;
              d->next = desc;
              desc = d;
              done_al = true;
            }
        // PI
        if (!done_pi) {
          if (strcmp(atm->resn, "HIS") == 0)
            if (strcmp(atm->tag, "CG") == 0 || strcmp(atm->tag, "ND1") == 0 || strcmp(atm->tag, "CD2") == 0 || 
                strcmp(atm->tag, "CE1") == 0 || strcmp(atm->tag, "NE2") == 0) 
              {
                d = new Descriptor(d_size++);
                d->crd = (atm->acid_start->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->crd) * (1.0/5.0);
                d->s->mnsp = 0x02;
                d->psurf = 0;
                //                d->colors.copy(atm->colors);
                d->atom = atm;
                d->next = desc;
                desc = d;
                done_pi = true;
              }
          if (strcmp(atm->resn, "PHE") == 0)
            if (strcmp(atm->tag, "CG") == 0 || strcmp(atm->tag, "CD1") == 0 || strcmp(atm->tag, "CD2") == 0 || 
                strcmp(atm->tag, "CE1") == 0 || strcmp(atm->tag, "CE2") == 0 || strcmp(atm->tag, "CZ") == 0) 
              {
                d = new Descriptor(d_size++);
                d->crd = (atm->acid_start->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->crd
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->next->crd) * (1.0/6.0);
                d->s->mnsp = 0x02;
                d->psurf = 0;
                //                d->colors.copy(atm->colors);
                d->atom = atm;
                d->next = desc;
                desc = d;
                done_pi = true;
              }
          if (strcmp(atm->resn, "TRP") == 0)
            if (strcmp(atm->tag, "CG") == 0 || strcmp(atm->tag, "CD1") == 0 || strcmp(atm->tag, "CD2") == 0 || 
                strcmp(atm->tag, "NE1") == 0 || strcmp(atm->tag, "CE2") == 0 || strcmp(atm->tag, "CE3") == 0 ||
                strcmp(atm->tag, "CZ2") == 0 || strcmp(atm->tag, "CZ3") == 0 || strcmp(atm->tag, "CH2") == 0) 
              {
                d = new Descriptor(d_size++);
                d->crd = (atm->acid_start->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->crd
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->next->crd
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->next->next->crd
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->next->next->next->crd
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->next->next->next->next->crd) 
                  * (1.0/9.0);
                d->s->mnsp = 0x02;
                d->psurf = 0;
                //                d->colors.copy(atm->colors);
                d->atom = atm;
                d->next = desc;
                desc = d;
                done_pi = true;
              }
          if (strcmp(atm->resn, "TYR") == 0)
            if (strcmp(atm->tag, "CG") == 0 || strcmp(atm->tag, "CD1") == 0 || strcmp(atm->tag, "CD2") == 0 || 
                strcmp(atm->tag, "CE1") == 0 || strcmp(atm->tag, "CE2") == 0 || strcmp(atm->tag, "CZ") == 0) 
              {
                d = new Descriptor(d_size++);
                d->crd = (atm->acid_start->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->crd 
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->crd
                          + atm->acid_start->next->next->next->next->next->next->next->next->next->next->crd) * (1.0/6.0);
                d->s->mnsp = 0x02;
                d->psurf = 0;
                //                d->colors.copy(atm->colors);
                d->atom = atm;
                d->next = desc;
                desc = d;
                done_pi = true;
              }
        }
      }
      atm = atm->next;
    }
  }
}

void Molecule::get_seq() {
  /* 
     Preberemo aminokislinsko zaporedje proteina (samo query verig, ki jih je lahko vec) in ga shranimo v sequence tipa map.
  */
  int resi = -XLARGE;
  char tmp_chain_id = ' ';
  Atom *at;

  for (at = atom; at != NULL; at = at->next) {
    if (chain_id.find(at->chain_id) != string::npos && !at->het && amino.find(at->resn) != string::npos) {
      if (at->resi != resi || at->chain_id != tmp_chain_id) {
        tmp_chain_id = at->chain_id;
        resi = at->resi;
        
        sequence.insert( new Residue(at->chain_id, one_letter_code(at->resn), at->resi) ); // sortirana po paru (chain_id, resi)
        
#ifdef VERB
        cout << "Molecule::get_seq() at->chain_id = " << at->chain_id << " at->resi = " << at->resi << " one_letter_code = " << one_letter_code(at->resn) << endl;
#endif
          
      }
    }
  }

  /* sedaj zapisemo minuse tam,kjer so vrzeli v zaporedju */ 
  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
    char chain_id1 = (*it)->chain_id;
    int resi1 = (*it)->resi;
    set<Residue*>::iterator it2 = it; it2++;
    if (it2 != sequence.end() ) {
      char chain_id2 = (*it2)->chain_id;
      int resi2 = (*it2)->resi;

      if (chain_id1 == chain_id2) 
        for (int i = resi1 + 1; i < resi2; i++) {
          sequence.insert ( new Residue(chain_id1, '-', i) );
#ifdef VERB
          cout << "Molecule::get_seq() missing residue chain_id1 = " << chain_id1 << " resi = " << i << endl;
#endif
        }
    }
  }

#ifdef VERB
  for (set<Residue*>::iterator it=sequence.begin(); it != sequence.end(); it++) 
    cout << "Molecule::get_seq() chain_id = " << (*it)->chain_id << " resn = " << (*it)->resn << " resi = " << (*it)->resi << endl;
#endif


}

void Molecule::binarize_cons_seq() {

  vector<int> bin;
  bin.clear();
  size_t max = 0;

  /* odpremo toliko bin-ov kolikor je najveckrat ohranjena aminokislina v zaporedju */
  /* vsak bin[i] predstavlja stevilo residujev, ki so ohranjeni i-krat */

  for (set<Residue*>::iterator it=sequence.begin(); it != sequence.end(); it++) {

#ifdef VERB
    int chain_id = it->chain_id; // chain_id
    int resi = it->resi; // resi
    Residue r(chain_id,resi);
    cout << "Molecule::binarize_cons_seq()  it->cons "<< it->cons << " sequence[" << chain_id << "," << resi << "]->cons = " << sequence.find(&r)->cons << endl;  
#endif

    if ((*it)->cons < max) { 
      bin[(*it)->cons]++;
    } 
    else {
      max = (*it)->cons;
      while (bin.size() < max + 1) bin.push_back(0);
      bin[max]++;
    }
    
  }

#ifdef VERB
  for (size_t j = 0; j < bin.size(); j++) {
    cout << "Molecule::binarize_cons_seq  bin["<<j<<"] = " << bin[j] << endl;
  }
#endif

  /* imamo 20 bin-ov, ceprav je koncna natancnost 10 */
  int num = 0, c = 19;
  float f = 0.05;

  for (size_t j = bin.size(); j > 0; j--) {
    num += bin[j-1];
#ifdef VERB
    cout << "Molecule::binarize_cons_seq num = "<<num<<endl;
    cout << "Molecule::binarize_cons_seq (float) num = "<<(float) num<<" (float) sequence.size() * f = "<< (float) sequence.size() * f <<endl;
    cout << "Molecule::binarize_cons_seq c = "<<c<<endl;
    cout << "Molecule::binarize_cons_seq f = "<<f<<endl;
#endif

    while ((float) num > (float) sequence.size() * f ) {
      c--;
      f += 0.05;
    }

    bin[j-1] = c ;

#ifdef VERB
    cout << "Molecule::binarize_cons_seq c = "<<c<<endl;
    cout << "Molecule::binarize_cons_seq f = "<<f<<endl;
    cout << "Molecule::binarize_cons_seq bin["<<j-1<<"] = " << bin[j-1] << endl;
    cout << endl;
#endif

  }

  for (set<Residue*>::iterator it=sequence.begin(); it != sequence.end(); it++) {
    (*it)->rawcons = (*it)->cons;  // shranimo nebinarizirane cons-e
    (*it)->cons = bin[(*it)->cons] / 2; // 20 -> 10
  }
}


void Molecule::out_info() {
  /*
    Izpisemo informacije o kakovosti alignmenta.
  */

  ofstream out_file((add_end_slash(OUTDIR) + "info.json").c_str());

  /* izpisemo score in drobnarije */
  out_file << "{";
//  out_file << "\"aligned_vertices\":" << SCONS << ",";
//  out_file << "\"fingerprint\":" << SIG_N << ",";
//  out_file << "\"surf_vector_angle\":" << SURF_VECTOR_ANGLE << ",";
//  out_file << "\"rmsd\":" << CALPHA_RMSD << ",";
//  out_file << "\"e_value\":" << scientific << setprecision(3) << BLOSUM_SCORE << ",";
  out_file << "\"z_score\":" << Z_SCORE << ",";
  out_file << "\"qpdb\":" << "\"" << pdb_id << "\"" << ",";
  out_file << "\"qcid\":" << "\"" << chain_id << "\"" << ",";

  /* izpisemo datoteke z bioloskimi assemblyi (prva je asymmetric unit) */
  out_file << "\"biof\":[\"" << pdb_id << ".cons.pdb" << "\"";
  
  for (map<pair<int,int>, BioUnit*>::iterator it = biounit.begin(); it != biounit.end(); it++) {
    if (it->second->filename.size() > 0) {
      out_file << ",";
      out_file << "\"" << it->second->filename << "\"";
    }
  }
  out_file << "]";
  out_file << "}";
  out_file << endl;

  out_file.close();
}



void Molecule::out_alignments(Subgraph *s, const string &json_file) {
  /*
    Izpisemo listo alignanih pdb-jev, njihovih dolzin in score-e za vsak posamezen klaster v JSON formatu
  */

  //~ ofstream out_file((add_end_slash(OUTDIR) + "alignments.json").c_str());
  ofstream out_file(json_file.c_str());
  /* v results file zapisemo json */
//  out_file<<"__REST:[";
  out_file<<"[";

  for (vector<Molecule*>::iterator it=lista_molekul.begin(); it != lista_molekul.end(); it++) {
	dbgmsg("pdb_id : " << (*it)->pdb_id);

    if (it != lista_molekul.begin())
      out_file<<",";
    
    /* izpisemo alignment score za vsak klaster */
    set<int> clus;  // ze obiskani klastri
    clus.clear();

    // dobimo stevilke klastrov, da bomo lahko sli po vrsti po klastrih
    for (multimap<Residue*,Residue*>::iterator it2 = (*it)->align.begin(); it2 != (*it)->align.end(); it2++) {
        clus.insert( it2->second->cd->cluster_id );
    }

    out_file<<"{";
    out_file<<"\"pdb_id\":\""<<(*it)->pdb_id<<"\",";
    out_file<<"\"chain_id\":\""<<(*it)->chain_id<<"\",";
    out_file<<"\"nfp\":"<<(*it)->nfp<<",";


    // ... pa se ime proteina (ki ga moramo prebrati iz ustrezne datoteke)
    string name = add_end_slash(LIGDIR) + ((*it)->pdb_id + (*it)->chain_id);
    //    string name = add_end_slash(LIGDIR) + "names/" + ((*it)->pdb_id + (*it)->chain_id);
    ifstream f(name.c_str());
    string text = "-";
    if (f.is_open()) { 
      getline (f, text);
    }
    f.close();
      
    /* escapamo posebne znake : " in \ */
    text = json_escape_special_characters(text);
    out_file << "\"name\":\"" << text <<"\",";


    //out_file<<"\"size\":"<<to_string((*it)->unique_seq_size)<<",";
    
    out_file<<"\"alignment\":[";

    // gremo po stevilkah klastrov
    for (set<int>::iterator itc = clus.begin(); itc != clus.end(); itc++) {
      if (itc != clus.begin()) {
	out_file<<",";
      }
      bool comma = false;
      /* gremo po ujemajocih se aminokislinah */
      for (multimap<Residue*,Residue*>::iterator it2 = (*it)->align.begin(); it2 != (*it)->align.end(); it2++) {

	Residue *r1 = it2->first;
	Residue *r2 = it2->second;
	// gremo po vrsti po klastrih (po narascajocih stevilkah)
	if (*itc == r2->cd->cluster_id) {
	  if (comma) {
	    out_file<<",";
	  }
	  else {
	    // izpisemo score za ta klaster
	    ClusterData *c = r2->cd;
	    out_file<<"{";  // en alignment
	    out_file<<"\"scores\":";
	    out_file<<"{";
	    out_file<<"\"alignment_no\":"<<c->cluster_id<<",";
	    out_file<<"\"aligned_vertices\":"<<c->scons<<",";
	    out_file<<"\"e_value\":"<<scientific<<setprecision(2)<<c->blosum_score<<",";
	    out_file<<"\"rmsd\":"<<fixed<<setprecision(1)<<c->calpha_rmsd<<",";
	    out_file<<"\"sva\":"<<fixed<<setprecision(2)<<(isnan(c->surf_vector_angle) ? 1.0e-2 : c->surf_vector_angle) <<",";
//	    out_file<<"\"sva\":"<<fixed<<setprecision(2)<<c->surf_vector_angle<<",";
	    out_file<<"\"z_score\":"<<fixed<<setprecision(2)<<z_score(c->cluster_score)<<",";
	    out_file<<"\"alignment_score\":"<<fixed<<setprecision(2)<<c->cluster_score;
	    out_file<<"}";
	    out_file<<",";

	    out_file<<"\"rotation_matrix\":";
	    out_file<<"[";
	    out_file<<"[";
	    out_file << (float) gsl_matrix_get(c->U, 0, 0) << ",";
	    out_file << (float) gsl_matrix_get(c->U, 0, 1) << ",";
	    out_file << (float) gsl_matrix_get(c->U, 0, 2);
	    out_file<<"],";
	    out_file<<"[";
	    out_file << (float) gsl_matrix_get(c->U, 1, 0) << ",";
	    out_file << (float) gsl_matrix_get(c->U, 1, 1) << ",";
	    out_file << (float) gsl_matrix_get(c->U, 1, 2);
	    out_file<<"],";
	    out_file<<"[";
	    out_file << (float) gsl_matrix_get(c->U, 2, 0) << ",";                 
	    out_file << (float) gsl_matrix_get(c->U, 2, 1) << ",";                 
	    out_file << (float) gsl_matrix_get(c->U, 2, 2);
	    out_file<<"]";
	    out_file<<"],";

	    out_file<<"\"translation_vector\":";
	    out_file<<"[";
	    out_file << (float) gsl_vector_get(c->t, 0) << ",";
	    out_file << (float) gsl_vector_get(c->t, 1) << ",";
	    out_file << (float) gsl_vector_get(c->t, 2);
	    out_file<<"],";

	    out_file<<"\"aligned_residues\":[";

	  }
	  comma = true;
//	  out_file << r1->resn << "." << r1->resi << "." << r1->chain_id;
//	  out_file << "-";
//	  out_file << r2->resn << "." << r2->resi << "." << r2->chain_id;
//	  out_file << r1->resi << ":" << r1->chain_id;
	  out_file<<"{";
	  out_file << "\"a\":\"" << r1->resn << ":" << r1->resi << ":" << r1->chain_id << "=" << r2->resn << ":" << r2->resi << ":" << r2->chain_id << "\",";
	  out_file << "\"c\":\"" << (r2->flx?"x":"") << (r1->sig_id?" f":"") << "\"";
	  out_file<<"}";
//	  out_file<<"{";
//	  out_file << "\"resn1\":\"" << r1->resn << "\",";
//	  out_file << "\"resi1\":" << r1->resi << ",";
//	  out_file << "\"chain_id1\":\"" << r1->chain_id << "\",";
////	  out_file << "\"flx\":" << (int) r2->flx << ",";
//	  out_file << "\"cl\":\"" << (r2->flx?"flx":"") << (r1->sig_id?" fp":"") << "\",";
//	  out_file << "\"resn2\":\"" << r2->resn << "\",";
//	  out_file << "\"resi2\":" << r2->resi << ",";
//	  out_file << "\"chain_id2\":\"" << r2->chain_id << "\"";
//	  out_file<<"}";
////	  out_file<<"\"";
//	  out_file << r1->resn << ":" << r1->resi << ":" << r1->chain_id;
//	  out_file << (r2->flx ? "-" : "=");
//	  out_file << r2->resn << ":" << r2->resi << ":" << r2->chain_id;
//	  out_file<<"\"";
	  
	}
      }
      out_file<<"]"; // ... aligned_residues
      out_file<<"}"; // .. konec "en alignment"
    }
    out_file<<"]"; // ... info
    out_file<<"}"; // ... alignan protein objekt

  }
  out_file<<"]";  // ... __REST
  out_file<<endl;
    
  out_file.close();

}







void Molecule::goodness_of_prediction() {
  /* 
     Izracuna specificity, sensitivity in significance of prediction.

     NOTE (JAN/26/2010) : Sedaj zadostuje ze, da je bfactor>BFACTOR enega samega atoma v aminokislin 
                          (prej so morali imeti vsi atomi bfactor razlicen od nic).
  */

  int tmp_resi = -9999999, tmp_bfactor_resi = -9999999;
  char tmp_chain = ' ';
  int Os = 0, Ps = 0, Ss = 0;
  int Is = 0;
  int Overlap = 0; // presek med dvema tipoma binding siteov, npr. med protein-protein in protein-small ligand
  map<char, int> Os_part;

  Os_part['x'] = 0;
  Os_part['y'] = 0;
  Os_part['z'] = 0;
  Os_part['w'] = 0;

  for (unsigned int i = 0; i < plain_pdb.size(); i++) {
    if (plain_pdb[i].compare(0, 4, "ATOM") == 0) {
      
      int resi = atoi( plain_pdb[i].substr(22,4).c_str() );

      float bfactor = atof( plain_pdb[i].substr(60,6).c_str() );

      char c_id = plain_pdb[i].at(21);

  
      if (chain_id.find(c_id) != string::npos) {


        /* this calculates the number of correctly predicted binding sites residues */
        if (tmp_resi != resi || tmp_chain != c_id) {
          tmp_chain = c_id;
          tmp_resi = resi;
          
          Ss++;
          
          /* velikost celokupnega interface-a se v vsakem primeru poveca */
          if ( imap.find(make_pair(c_id, resi)) != imap.end() ) Is++;

          /* koliko residujev je v preseku med dvema tipoma b.s.-jev ? */
          if ( (int) imap.count(make_pair(c_id, resi)) > 1 ) Overlap++;
          
          
        }



        /* velikost correctly predicted (Os) izracunamo samo v primeru, da ima ak ugoden bfactor */
        
        if (bfactor > BFACTOR) {
          
          if (tmp_bfactor_resi != resi) {
            tmp_bfactor_resi = resi;
            
            Ps++;
            if ( imap.find(make_pair(c_id, resi)) != imap.end() ) {

              Os++;

              /* ker lahko vsaka query residue pripada raznim binding sitesom */
              pair<multimap<pair<char,int>,char>::iterator, multimap<pair<char,int>,char>::iterator> ret = imap.equal_range(make_pair(c_id, resi));
              /* moramo po vseh razlicnih b.s.-jih katerim en residue lahko pripada */
              for (multimap<pair<char,int>,char>::iterator it = ret.first; it != ret.second; it++) {
//              char veriga = imap.find(make_pair(c_id, resi))->second;
                char veriga = it->second;
//              Os_part[veriga]++;
                if (veriga == '1')
                  Os_part['w']++; // small ligand
//                else if (isalpha(veriga))
                else if (veriga == '@')
                  Os_part['x']++;  // protein-protein
                else if (veriga == '2')
                  Os_part['y']++;  // dna
                else if (veriga == '*')
                  Os_part['z']++;  // fn import
              }

            }
          }
        }

        
      }
    }
  }
  
  //  int Is = 0;
  map<char, int> Is_part;

  Is_part['x'] = 0;
  Is_part['y'] = 0;
  Is_part['z'] = 0;
  Is_part['w'] = 0;

  /* dobimo velikost interfacea samo na query chainu */
  for (multimap<pair<char,int>, char>::iterator it = imap.begin(); it != imap.end(); it++) {
    if (chain_id.find(it->first.first) != string::npos) { 
//      Is++;
      char veriga = it->second;
//      Is_part[it->second]++;

      if (veriga == '1')
        Is_part['w']++; // small ligand
//      else if (isalpha(veriga))
      else if (veriga == '@')
        Is_part['x']++;  // protein-protein
      else if (veriga == '2')
        Is_part['y']++;  // dna
      else if (veriga == '*')
        Is_part['z']++;  // fn import

//      /* presek med small-ligand in protein-protein */
//      if (veriga == '1') {
//        int resi = it->first->second;
//        presek[resi] = 1;
//      }

    }
  }

  cout << "TITL>\tPROTEIN_ID\tTOTAL\tINTER\tSIGNIFI\tPREDI\tSPECI\tSENSI\tCORRE";

  /* izpisemo tudi posamicne sensitivity za vsako verigo posebej */
  for (map<char,int>::iterator it = Os_part.begin(); it != Os_part.end(); it++) {
    cout << "\t" << it->first;
  }
  /* in se stevilo residujev, ki pripadajo dvema ali vec b.s.-teom */
  cout << "\tOVERL" << endl;

  printf("DATA>\t%s\t%d\t%d\t%.1e\t%d\t%.3f\t%.3f\t%d", pdb_id.c_str(), Ss, Is, significance(Ss, Ps, Is, Os), Ps, (float) Os / Ps, (float) Os / Is, Os);
  
  for (map<char,int>::iterator it = Os_part.begin(); it != Os_part.end(); it++) {
    printf("\t%.3f", (float) it->second / Is_part[it->first]);
//    printf("\t%.3f", (float) Is_part[it->first]);
  }

  printf("\t%d", Overlap);

//  for (map<char,int>::iterator it = Os_part.begin(); it != Os_part.end(); it++) {
//    printf("\t%.3f", (float) Is_part[it->first]);
//  }

  cout << endl;
}






//	void Molecule::conservation0(float percent) {
//	  /*
//	    Dolocimo ohranjenost vseh aminokislin v query zaporedju.
//	  */
//	  
//	  /* gremo po zaporedju mol1 - query proteina */
//	  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
//	
//	    /* za vsako mesto v query zaporedju inicializiramo nov graf v item in qmax */
//	    if ((*it)->resn != '-') {
//	
//	#ifdef VERB
//	      cout << "Molecule::conservation Sequence position " << (*it)->resn << (*it)->resi << (*it)->chain_id << endl;
//	#endif 
//	
//	      /* postavimo na nic, ker proceduro veckrat pozenemo */
//	      (*it)->cons = 0;
//	      /* gremo po vseh prileganih zaporedjih */
//	      for (vector<Molecule*>::iterator it2 = lista_molekul.begin(); it2 != lista_molekul.end(); it2++) {
//	        
//	#ifdef VERB
//	        cout << "Molecule::conservation "<< (*it2)->unique_seq_size << " " << sequence.size() << endl;
//	#endif
//	
//	        /* upostevamo le prva zaporedja, ki imajo alignment length vecji od N% query zaporedja */
//	        if ( (float) (*it2)->unique_seq_size > (float) (percent * sequence.size()) ) { 
//	
//	          int chain_id = (*it)->chain_id;
//	          int resi = (*it)->resi;
//	#ifdef VERB
//	          int main_id = (*it2)->main_id;
//	          cout << "Molecule::conservation chain_id = " << chain_id << " resi = " << resi << " main_id = " << main_id << endl;
//	#endif
//	
//	          Residue r(chain_id,resi);
//	          pair<multimap<Residue*, Residue*>::iterator, multimap<Residue*, Residue*>::iterator> ret = (*it2)->align.equal_range(&r);
//	
//	          /* ce je i-to mesto ohranjeno v vsaj enem klastru, povecamo stevec ohranjenosti mesta i v query zaporedju */
//	          if (ret.first != ret.second) {
//	            (*it)->cons++; 
//	            
//	          }
//	
//	        }
//	      }
//	    }
//	  }
//	  
//	}




//void Molecule::conservation2(float cluster_score) {
//  /*
//    Dolocimo ohranjenost vseh aminokislin v query zaporedju. Upostevamo samo zaporedja, ki imajo 
//    first_cluster_score vecji od neke vrednosti.
//  */
//  
//  /* gremo po zaporedju mol1 - query proteina */
//  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
//    
//    /* za vsako mesto v query zaporedju inicializiramo nov graf v item in qmax */
//    if ((*it)->resn != '-') {
//      
//#ifdef VERB
//      cout << "Molecule::conservation Sequence position " << (*it)->resn << (*it)->resi << (*it)->chain_id << endl;
//#endif 
//      
//      /* postavimo na nic, ker proceduro veckrat pozenemo */
//      (*it)->cons = 0;
//      /* gremo po vseh prileganih zaporedjih */
//      for (vector<Molecule*>::iterator it2 = lista_molekul.begin(); it2 != lista_molekul.end(); it2++) {
//        
//#ifdef VERB
//        cout << "Molecule::conservation "<< (*it2)->first_cluster_score << " "<< (*it2)->unique_seq_size << " " << sequence.size() << endl;
//#endif                                                                                                     
//        
//        /* upostevamo le prva zaporedja, ki imajo first_cluster_score vecji od neke mejne vrednosti */
//        if ((*it2)->first_cluster_score > cluster_score) { 
//          
//          int chain_id = (*it)->chain_id;
//          int resi = (*it)->resi;
//#ifdef VERB
//          int main_id = (*it2)->main_id;
//          cout << "Molecule::conservation chain_id = " << chain_id << " resi = " << resi << " main_id = " << main_id << endl;
//#endif
//          
//          Residue r(chain_id,resi);
//          pair<multimap<Residue*, Residue*>::iterator, multimap<Residue*, Residue*>::iterator> ret = (*it2)->align.equal_range(&r);
//          
//          /* ce je i-to mesto ohranjeno v vsaj enem klastru, povecamo stevec ohranjenosti mesta i v query zaporedju */
//          if (ret.first != ret.second) {
//            (*it)->cons++; 
//          }
//        }
//      }
//    }
//  }
//}


void Molecule::conservation(float z) {
  /*
    Dolocimo ohranjenost vseh aminokislin v query zaporedju. Upostevamo samo zaporedja, ki imajo 
    iz first_cluster_score-a izracunani z_score > z.
  */
  
  /* gremo po zaporedju mol1 - query proteina */
  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
    
    /* za vsako mesto v query zaporedju inicializiramo nov graf v item in qmax */
    if ((*it)->resn != '-') {
      
#ifdef VERB
      cout << "Molecule::conservation Sequence position " << (*it)->resn << (*it)->resi << (*it)->chain_id << endl;
#endif 
      
      /* postavimo na nic, ker proceduro veckrat pozenemo */
      (*it)->cons = 0;
      /* gremo po vseh prileganih zaporedjih */
      for (vector<Molecule*>::iterator it2 = lista_molekul.begin(); it2 != lista_molekul.end(); it2++) {
        
#ifdef VERB
        cout << "Molecule::conservation "<< (*it2)->first_cluster_score << " "<< (*it2)->unique_seq_size << " " << sequence.size() << endl;
#endif                                                                                                     
        
        /* upostevamo le zaporedja, ki imajo z_score za first_cluster_score vecji od neke mejne vrednosti (z) */
        if (z_score((*it2)->first_cluster_score) > z) { 
          
          int chain_id = (*it)->chain_id;
          int resi = (*it)->resi;
#ifdef VERB
          int main_id = (*it2)->main_id;
          cout << "Molecule::conservation chain_id = " << chain_id << " resi = " << resi << " main_id = " << main_id << endl;
#endif
          
          Residue r(chain_id,resi);
          pair<multimap<Residue*, Residue*>::iterator, multimap<Residue*, Residue*>::iterator> ret = (*it2)->align.equal_range(&r);
          
          /* ce je i-to mesto ohranjeno v vsaj enem klastru, povecamo stevec ohranjenosti mesta i v query zaporedju */
          if (ret.first != ret.second) {
            (*it)->cons++; 
          }
        }
      }
    }
  }
}

void Molecule::fingerprint() {
  /*
    Dolocimo fingerprint aminokisline.
  */

  int max_cons = 8;

  vector<pair<ChResi,int> > sig_res; // signature residues


  /* gremo po zaporedju mol1 - query proteina */
  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {

    /* dolocimo najbolj ohranjene */
    if ((*it)->resn != '-' && (signed) (*it)->cons > max_cons) {
      max_cons = (*it)->cons;
    }

  }

  /* gremo po zaporedju mol1 - query proteina */
  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
    /* potencialne signature ak so kar vse, ki so zelo ohranjene (binned cons >= max_cons) */ 
    if ((*it)->resn != '-' && (signed) (*it)->cons > max_cons - 1) {
#ifdef VERB
      cout << "signature_residues2 " << (*it)->chain_id << " " << (*it)->resi << " "  << (*it)->cons << " " << endl;
#endif
      sig_res.push_back(make_pair( ChResi((*it)->chain_id, (*it)->resi), (*it)->cons));
    }
  }

  if (sig_res.empty() ) {  // ce ni nobenega sig_res, potem pustimo samo 100 prvih prileganih zaporedij
    for (unsigned int i = 100; i < lista_molekul.size(); i++) {
      delete lista_molekul[i];
      lista_molekul.erase(lista_molekul.begin() + i );
      i--;
    }
  }
  else {  // ce najdemo sig res, potem izbrisemo prilegana zaporedja, ki jih nimajo dovolj

    /* v query zaporedju oznacimo signature residues */
    for (unsigned int i = 0; i < sig_res.size(); i++) {
      char chain_id0 = sig_res[i].first.first;
      int resi0 = sig_res[i].first.second;
      Residue r(chain_id0, resi0);
      (*sequence.find(&r))->sig_id = 1; // 1 oznacuje, da je to fingerprint residue
    }

    /* dolocimo stevilo signature residuejev v vsakem prileganem zaporedju */
    for (unsigned int i = 0; i < lista_molekul.size(); i++) {
      int num_sig = 0;
      /* preverimo ce to zaporedje ima SIG_N od vseh kolikor jih upostevamo (vse v SIG_RESI_CHECKED) signature residues */
      for (unsigned int j = 0; j < sig_res.size(); j++) {

        char chain_id0 = sig_res[j].first.first;
        int resi0 = sig_res[j].first.second;
        Residue r(chain_id0, resi0);
        pair<multimap<Residue*,Residue*>::iterator, multimap<Residue*,Residue*>::iterator> ret = lista_molekul[i]->align.equal_range(&r);
        
	bool has_sig = false;
        /* ce je to signature mesto ohranjeno, povecamo sig_stevec tega zaporedja in oznacimo v zaporedju */
	if (ret.first != ret.second) {
	  for (multimap<Residue*,Residue*>::iterator it = ret.first; it != ret.second; it++) {
	    if (!it->second->flx) has_sig = true;  // za fingerprint residueje priznavamo samo non-flexible residueje
          }
        }
	if (has_sig) num_sig++;

      }
      // stevilo fingerprint residuejev za ta alignment
      lista_molekul[i]->nfp = num_sig;
#ifdef VERB
      cout << "Molecule::fingerprint() nfp=" << num_sig << endl;
#endif
      /* ce ni vseh SIG_N (=NFP) signature residues, potem izbrises to zaporedje */
      if ( num_sig < (int) ( sig_res.size()*NFP ) && lista_molekul.size() > 1) {
	delete lista_molekul[i];
	lista_molekul.erase(lista_molekul.begin() + i );
	i--;
      }
    }
  }
}



void Molecule::statistics() {
  /* gremo po zaporedju mol1 - query proteina */
  for (set<Residue*>::iterator it = sequence.begin(); it != sequence.end(); it++) {
    /* za vsako mesto v query zaporedju inicializiramo nov graf v item in qmax */
    if ((*it)->resn != '-') {
      ofstream statf ((to_string((*it)->resi) + ".stat").c_str());
      
#ifdef VERB
      cout << "Molecule::statistics Sequence position " << (*it)->resn << (*it)->resi << (*it)->chain_id << endl;
#endif 
      /* gremo po vseh prileganih zaporedjih */
      for (vector<Molecule*>::iterator it2 = lista_molekul.begin(); it2 != lista_molekul.end(); it2++) {
        
#ifdef VERB
        cout << "Molecule::statistics "<< (*it2)->first_cluster_score << " "<< (*it2)->unique_seq_size << " " << sequence.size() << endl;
#endif                                                                                                     
        
//        /* upostevamo le zaporedja, ki imajo z_score za first_cluster_score vecji od neke mejne vrednosti (z) */
//        if (z_score((*it2)->first_cluster_score) > z) { 
          
	int chain_id = (*it)->chain_id;
	int resi = (*it)->resi;
#ifdef VERB
	int main_id = (*it2)->main_id;
	cout << "Molecule::statistics chain_id = " << chain_id << " resi = " << resi << " main_id = " << main_id << endl;
#endif
	
	Residue r(chain_id,resi);
	pair<multimap<Residue*, Residue*>::iterator, multimap<Residue*, Residue*>::iterator> ret = (*it2)->align.equal_range(&r);
	if (ret.first != ret.second) {
	  Residue *rmax = max_element(ret.first, ret.second, by_noflx_cluster_score)->second;
	  statf << rmax->cd->cluster_score << endl;
	}
      }
      statf.close();
    }
  }
}
