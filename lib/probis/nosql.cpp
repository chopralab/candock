#include "nosql.h"
#include "molecule.h"
#include "residue.h"
#include "clusterdata.h"
#include "denorm.h"
#include "item.h"
#include "subgraph.h"
#include "product.h"
#include "desc.h"
#include "atom.h"
#include "debug.hpp"

#ifndef _WINDOWS

#include <sys/file.h> /* for flock(2) */
#include <sys/stat.h> /* for S_* constants */
#include <unistd.h> /* for close(2) prototypes */

#endif

#include <string.h> /* for strerror(3) prototype */
#include <stdio.h> /* for fprintf(3),printf(3),stderr protype */
#include <errno.h> /* for errno prototype */


//#include "kabsch.h"

// sortiramo po cluster_score prvega klastra (first_cluster_score je cluster_score prvega klastra, najvecji cluster_score), nato, ce imamo enakost, po dolzini celotnega alignmenta
bool by_cluster_score_unique_seq_length(Molecule *i, Molecule *j) {
  return i->first_cluster_score > j->first_cluster_score ? true : (i->first_cluster_score == j->first_cluster_score && i->unique_seq_size > j->unique_seq_size ? true : false); 
}

//~ NoSql::NoSql(Molecule *m) {
  //~ this->mol1 = m;
  //~ if (m != 0) set_db_file(m);
//~ //  this->db_file = NOSQL_DIR + "/" + m->pdb_id + m->chain_id + ".nosql";
//~ }
//~ 
//~ NoSql::NoSql(const string& mol_pdb_id, const string& mol_chain_id) {
  //~ this->mol1 = 0;
  //~ set_db_file(mol_pdb_id, mol_chain_id);
//~ }

void NoSql::set_db_file(Molecule *m) {
//  this->db_file = NOSQL_DIR + "/" + m->pdb_id + m->chain_id + ".nosql";
  this->db_file = m->pdb_id + m->chain_id + ".nosql";
}

void NoSql::set_db_file(const string& mol_pdb_id, const string& mol_chain_id) {
//	db_file = NOSQL_DIR + "/" + mol_pdb_id + mol_chain_id + ".nosql";
	db_file = mol_pdb_id + mol_chain_id + ".nosql";
}


void NoSql::append_file (string row) {
  /*
    Appendamo vrstico na konec baze.
  */
//TODO: Fix later
#ifndef _WINDOWS

  int fd;
  fd = open((add_end_slash(OUTDIR) + db_file).c_str(),O_RDWR|O_CREAT,S_IRUSR|S_IWUSR); 

  /* Aquire an exclusive lock */ 
  if (flock(fd,LOCK_EX) == -1) {
  } 

#endif

  ofstream fn ((add_end_slash(OUTDIR) + db_file).c_str() , ios_base::app);
  if (fn.is_open()) {
    fn << row;
    fn.close();

#ifndef _WINDOWS

    if (flock(fd,LOCK_UN)==-1) {
    }
    if (close(fd)==-1) {
    }

#endif


  }
  else {
    throw Err("Error (NOSQL) : Unable to open file for appending", 17); 
  }
}

void NoSql::read_file (vector<string> &rows) {
  /*
    Preberemo cel .nosql file.
  */
  string line;
  ifstream fn ((add_end_slash(INDIR) + db_file).c_str());
  if (fn.is_open()) {
    while ( fn.good() ) {
      getline (fn,line);
      rows.push_back(line);
    }
    fn.close();
  }
  else {
    throw Err("Error (NOSQL) : Unable to open file", 18); 
  }
}

void NoSql::read_file (Molecule *m2, string &row) {
  /*
    Preberemo samo eno vrstico v .nosql datoteki.
  */
  string line;
  Denorm *denorm = new Denorm();
  ifstream fn ((add_end_slash(INDIR) + db_file).c_str());
  if (fn.is_open()) {
    bool ex = false;
    while ( fn.good() && !ex) {
      getline (fn,line);
      if (denorm->read_one(line)) {
	if (m2->pdb_id.compare(denorm->get_pdb_id()) == 0 && m2->chain_id.compare(denorm->get_chain_id()) == 0) {
	  row = line;
	  ex = true;
	}
      }
    }
    fn.close();
  }
  else {
    throw Err("Error (NOSQL) : Unable to open file", 19); 
  }
  delete denorm;
}

void NoSql::read() {
  /*
    Beremo iz nosql baze, vsak non-redundant protein je svoja datoteka v tej bazi. Inicializiramo vse
    spremenljivke v classu Molecule.
  */
  vector<string> rows;
  read_file(rows);

  for (vector<string>::iterator it = rows.begin(); it != rows.end(); it++) {
    Denorm *denorm = new Denorm();
    // ce je prebrana vrstica v pravilnem formatu in ce je sploh kaksen dober alignment ...
	dbgmsg("before if " << *it);
    if (denorm->read(*it, !_pairwise)) {
		dbgmsg(*it);
      this->mol1->lista_molekul.push_back(new Molecule(denorm->get_pdb_id() + ".pdb", denorm->get_chain_id()) );
      // preberemo vse klastre ...
      denorm->reset_cluster();
      ClusterData *c = NULL;
      //        cout << denorm->get_pdb_id() << "c=" << (c == NULL ? "null" : "not") << endl;
      while ((c = denorm->get_next_cluster()) != NULL) {
	// ... naknadno dolocimo pointer na Molecule*, ki pove, kateremu alignanemu proteinu pripada ta klaster
	c->m = this->mol1->lista_molekul.back();
	
	Residue *r1 = NULL, *r2 = NULL;
	denorm->reset_residue(c->cluster_id); // postavimo se na prvi par alignanih ak tega klastra
	
	// preberemo vse alignane residue-je tega klastra ...
	while (denorm->get_next_residue(r1, r2)) {
	  set<Residue*>::iterator its = this->mol1->sequence.find(r1);  // sortirane so po chain_id , nato resi, nato m (m tukaj nima veze)
	  if (its == this->mol1->sequence.end()) {
	    throw Err("Error (NOSQL)  : read()  : Cannot find query residue chain_id=" + to_string(r1->chain_id) + " resi=" + to_string(r1->resi) + ".", 21);
	  }
	  delete r1;
	  // ... tudi tu moramo naknadno dolociti pointer na Molecule*
	  r2->m = this->mol1->lista_molekul.back();
	  // ta f-ja mora biti nujno za read_sequence
	  this->mol1->lista_molekul.back()->align.insert(make_pair(*its, r2));
	}
      }
      
      /* cluster_score prvega klastra (ni nujno 0-ti !!!) zapisemo v Molecule::first_cluster_score */
      denorm->reset_cluster();
      this->mol1->lista_molekul.back()->first_cluster_score = denorm->get_next_cluster()->cluster_score;
      /* unique_seq_size je dolzina unique prileganega zaporedja (brez prekrivajocih se delov klastrov) */
      this->mol1->lista_molekul.back()->unique_seq_size = denorm->get_unique_seq_size();
      
    }
    delete denorm;
  }
  
  /* po first_cluster_score, nato znotraj enakega cluster_score-a po unique_seq_length */
  sort(this->mol1->lista_molekul.begin(), this->mol1->lista_molekul.end(), by_cluster_score_unique_seq_length);
	for (vector<Molecule*>::iterator it=this->mol1->lista_molekul.begin(); it != this->mol1->lista_molekul.end(); it++) {
		dbgmsg("pdb_id : " << (*it)->pdb_id);
	}

}

void NoSql::get_rota_trans_resi(Molecule *mol2, Item *one_clq) {
  /*
    Preberemo rotacijsko matriko in translacijski vektor za ALIGNMENT_NO-ti klaster prileganih proteinov
    v mol1 in mol2. Ce je reverse je true, potem je to oznaka, da sta mol1 in mol2 zamenjana.
  */
  cout << "NOSQL> Getting the rotation-translation matrix and member residues of cluster.." << endl;  
  string row;
  Denorm *denorm = new Denorm();

  try {
    // preberemo samo eno vrstico v .nosql datoteki, tisto, kjer je mol2->pdb_id in mol2->chain_id
    read_file(mol2, row);
    // ce je prebrana vrstica v pravilnem formatu ...
    if (!denorm->read(row, false)) throw Err("Error (NOSQL) : Cannot read .nosql file.", 1);
    // preberemo rotacijsko matriko in translacijski vektor
    one_clq->U = denorm->get_U(ALIGNMENT_NO);
    if (one_clq->U == NULL) throw Err("Error (NOSQL) : Cluster number --alno " + to_string(ALIGNMENT_NO) + " is out of range.", 1);
    one_clq->t = denorm->get_t(ALIGNMENT_NO);
    // dobimo conserved aminokisline na enem in drugem alignanem proteinu
    one_clq->cons = denorm->get_cons(ALIGNMENT_NO);
  } catch (Err e) {
    denorm->free(ALIGNMENT_NO);
    delete denorm;
    throw e;
  }
  denorm->free(ALIGNMENT_NO);
  delete denorm;
  
}


void NoSql::write(Molecule *mol2, Subgraph *subgraph, Score *score) {
	/*
    Izpisemo rezultate v NoSQL bazo...
  */
	NoSql::WriteStruct temp;
	writeStruct(mol2, subgraph, score, temp);
	append_file(temp.getString(2));
    
  // ce delamo ProBiS-Database, potem moramo izpisati se inverzno datoteko...
  if (_database) {
    set_db_file(mol2);
    append_file(temp.getString(1));
  }
}


void NoSql::write(const NoSql::WriteStruct& ws) {
	/*
    Izpisemo rezultate v NoSQL bazo...
  */
	append_file(ws.getString(2));
    
  // ce delamo ProBiS-Database, potem moramo izpisati se inverzno datoteko...
  if (_database) {
    set_db_file(ws.mol2_pdb_id, ws.mol2_chain_id);
    append_file(ws.getString(1));
  }
}


void NoSql::writeStruct(Molecule *mol2, Subgraph *subgraph, Score *score, NoSql::WriteStruct& tempStruct) {
  /*
    Izpisemo rezultate v zacasno strukturo, ki jo kasneje zapisemo v nosql bazo...
  */

  // izpisujemo samo, ce imamo kaj ...
  if (subgraph->clus.size() > 0) {
//    stringstream prefix;
//    prefix.str("");
    stringstream suffix;
    suffix.str("");

    // izpisemo ime alignanega proteina
//    prefix << mol2->pdb_id << "\t" << mol2->chain_id << "\t";
    
    // dolocimo stevilo klastrov za izpis
    size_t size;
    if (subgraph->clus.size() > (unsigned int) CLUS_TO_OUTPUT ) size = CLUS_TO_OUTPUT; else size = subgraph->clus.size(); // output only first N clusters
    
    // izpisemo cluster
    for (size_t i = 0; i < size; i++) {
      suffix << (i == 0 ? "(" : ",(");
      suffix << i; // cluster_id
      suffix << "," << subgraph->clus[i]->vert.size();
      suffix << "," << subgraph->clus[i]->N;  // NUMBER OF CLUSTERED CLIQUES
      suffix << "," << subgraph->clus[i]->calpha_rmsd;
      suffix << "," << subgraph->clus[i]->score_probe;
      suffix << "," << subgraph->clus[i]->score_descriptor_ratio;
      suffix << "," << subgraph->clus[i]->surf_vector_angle;
      suffix << "," << subgraph->clus[i]->blosum_score;
      suffix << "," << subgraph->clus[i]->cluster_score;
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 0, 0);
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 0, 1);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 0, 2);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 1, 0);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 1, 1);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 1, 2);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 2, 0);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 2, 1);                 
      suffix << "," << (float) gsl_matrix_get(subgraph->clus[i]->U, 2, 2);                 
      suffix << "," << (float) gsl_vector_get(subgraph->clus[i]->t, 0);                 
      suffix << "," << (float) gsl_vector_get(subgraph->clus[i]->t, 1);                 
      suffix << "," << (float) gsl_vector_get(subgraph->clus[i]->t, 2);                 
      suffix << ")";
    }
    
    suffix << "\t";
    
    // izpisemo result uniquely (ce sta enaka flx in ne-flx, potem izpisemo samo ne-flx), zanasamo se na to, da so ne-flexibini pred flexibilnimi
    for (size_t i = 0; i < size; i++) {
      suffix << (i == 0 ? "" : ",");
      set<string> unique;
      unique.clear();
      for (size_t j = 0; j < subgraph->clus[i]->vert.size(); j++) {

	stringstream s;
	s.str("");
//        s << (subgraph->clus[i]->vert[j]->row->desc2->num == -1 ? 'Y' : 'N');
        s << "," << subgraph->clus[i]->vert[j]->row->desc1->atom->chain_id;
        s << "," << one_letter_code(subgraph->clus[i]->vert[j]->row->desc1->atom->resn);
        s << "," << subgraph->clus[i]->vert[j]->row->desc1->atom->resi;
        s << "," << subgraph->clus[i]->vert[j]->row->desc2->atom->chain_id;
        s << "," << one_letter_code(subgraph->clus[i]->vert[j]->row->desc2->atom->resn);
        s << "," << subgraph->clus[i]->vert[j]->row->desc2->atom->resi;
        s << "," << i;
	// preverimo unikatnost kombinacije alignanih residuejev
	if (unique.insert(s.str()).second == true) {
	  suffix << (unique.size() == 1 ? "(" : ",(");
	  suffix << (subgraph->clus[i]->vert[j]->row->desc2->num == -1 ? 'Y' : 'N') << s.str(); // ce je ze bil dodan par alignanih ak z N, potem ne bo se enkrat dodan isti par z Y
	  suffix << ")";
	}

      }
    }
    suffix << "\n";
		tempStruct.suffix = suffix.str();
  }
	tempStruct.mol2_pdb_id = mol2->pdb_id; tempStruct.mol2_chain_id = mol2->chain_id;
	tempStruct.mol1_pdb_id = mol1->pdb_id; tempStruct.mol1_chain_id = mol1->chain_id;
}
