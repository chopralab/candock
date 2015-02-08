#include "denorm.h"
#include "geo.h"
#include "kabsch.h"
#include "clusterdata.h"
#include "residue.h"

gsl_matrix* Denorm::get_U(int cluster_id) { 
//  gsl_matrix *ret;
//  gsl_matrix_memcpy(ret, cluster[cluster_id]->U);
//  return ret; 
  ClusterData *c = cluster[cluster_id];
  return (c == NULL ? NULL : cluster[cluster_id]->U); 
//  return cluster[cluster_id]->U; 
}

gsl_vector* Denorm::get_t(int cluster_id) { 
//  gsl_vector *ret;
//  gsl_vector_memcpy(ret, cluster[cluster_id]->t);
//  return ret; 
  return cluster[cluster_id]->t; 
}

vector<pair<ChResi, int> > Denorm::get_cons(int cluster_id) { 
  /*
    Pretvorimo aligned za en klaster v drugacen format.
  */
  vector<pair<ChResi, int> > ret;
  reset_residue(cluster_id);
  Residue *r1 = NULL, *r2 = NULL;
  while(get_next_residue(r1, r2)) {
    ret.push_back(make_pair(ChResi(r1->chain_id, r1->resi), 1) );  // 1 je oznaka za prvi protein
    ret.push_back(make_pair(ChResi(r2->chain_id, r2->resi), 2) );  // 2 je oznaka za drugi protein
  }
  return ret;
}

void Denorm::free(int cluster_id) {
  cout << "Deleting denorm internal structures, except for U and t of cluster_id = " << cluster_id << " ..." << endl;
  for (map<int, ClusterData*>::iterator it = cluster.begin(); it != cluster.end(); it++) {
    if (it->first != cluster_id)
      delete it->second;
  }
  for (multimap<int, pair<Residue*, Residue*> >::iterator it = aligned.begin(); it != aligned.end(); it++) {
    delete it->second.first;
    delete it->second.second;
  }
}

bool Denorm::read_one(string &row) {
  /*
    Preberemo samo pdb_id in chain_id v denormalizirani vrstici in podatke zapisemo v interne strukture.
  */
  vector<string> all = reg_split(row, "\t");
  if (all.size() != 4) {
    cout << "Warning (DENORM) : Wrong format of this nosql line:" << row << endl;
    return false;
  }
  
  string pdb_id2 = all[0];
  string chain_id2 = (all[1].at(0) == '*' ? all[1].substr(1) : all[1]); // ce je inverted
  this->pdb_id = pdb_id2;
  this->chain_id = chain_id2;
  return true;
}

bool Denorm::read(string &row, bool filter) {
  /*
    Preberemo denormalizirano vrstico in podatke zapisemo v interne strukture. Hkrati tudi filtriramo 
    ven klastre s slabimi score-i.
  */
  vector<string> all = reg_split(row, "\t");
  if (all.size() != 4) {
    cout << "Warning (DENORM) : Wrong format of this nosql line:" << row << endl;
    return false;
  }
  
  string pdb_id2 = all[0];
  bool swap = (all[1].at(0) == '*' ? true : false);   // ce je zvezdica, sta protein1 in protein2 zamenjana 
  string chain_id2 = (all[1].at(0) == '*' ? all[1].substr(1) : all[1]); // ce je inverted
  string clusters = all[2];
  string results = all[3];


  this->pdb_id = pdb_id2;
  this->chain_id = chain_id2;

  // preberemo clusters
  vector<string> all_clusters = reg_split(clusters, "),(", "(", ")");
  for (vector<string>::iterator it = all_clusters.begin(); it != all_clusters.end(); it++) {
    vector<string> cluster_parts = reg_split(*it, ",");
    int cluster_id = atoi(cluster_parts[0].c_str());
    int scons = atoi(cluster_parts[1].c_str());
    float calpha_rmsd = atof(cluster_parts[3].c_str());
    float surf_vector_angle = atof(cluster_parts[6].c_str());
    double blosum_score = atof(cluster_parts[7].c_str());
    float cluster_score = atof(cluster_parts[8].c_str());

//    if (!filter || (blosum_score < BLOSUM_SCORE && calpha_rmsd < CALPHA_RMSD && surf_vector_angle < SURF_VECTOR_ANGLE && scons > SCONS)) {
    if (!filter || (z_score(cluster_score) > Z_SCORE)) {
      gsl_matrix *U = gsl_matrix_alloc(3, 3);
      gsl_vector *t = gsl_vector_alloc(3);
      
      gsl_matrix_set(U, 0, 0, atof(cluster_parts[9].c_str()));
      gsl_matrix_set(U, 0, 1, atof(cluster_parts[10].c_str()));
      gsl_matrix_set(U, 0, 2, atof(cluster_parts[11].c_str()));
      gsl_matrix_set(U, 1, 0, atof(cluster_parts[12].c_str()));
      gsl_matrix_set(U, 1, 1, atof(cluster_parts[13].c_str()));
      gsl_matrix_set(U, 1, 2, atof(cluster_parts[14].c_str()));
      gsl_matrix_set(U, 2, 0, atof(cluster_parts[15].c_str()));
      gsl_matrix_set(U, 2, 1, atof(cluster_parts[16].c_str()));
      gsl_matrix_set(U, 2, 2, atof(cluster_parts[17].c_str()));
      
      gsl_vector_set(t, 0, atof(cluster_parts[18].c_str()));
      gsl_vector_set(t, 1, atof(cluster_parts[19].c_str()));
      gsl_vector_set(t, 2, atof(cluster_parts[20].c_str()));
      
      // ... ce sta proteina zamenjana, potem invertiramo matriko in negiramo vektor
      if (swap) {
	Kabsch k;
	k.invert(U, t);
      }

      // pointer na alignan protein (Molecule*) dodamo kasneje (glej NoSql::read)
      cluster[cluster_id] = new ClusterData(cluster_id, scons, blosum_score, calpha_rmsd, surf_vector_angle, cluster_score, NULL, U, t);
    }

  }

  // ce ni noben klaster dovolj dober ...
  if (cluster.empty()) {
    return false;
  }

  // preberemo results
  set<ChResi> rid;
  // drugacen DB query, zato je vedno treba preverjati, po katerem proteinu smo dejansko iskali
  vector<string> all_results = reg_split(results, "),(", "(", ")");
  for (vector<string>::iterator it = all_results.begin(); it != all_results.end(); it++) {
      vector<string> result_parts = reg_split(*it, ",");
      char flx = result_parts[0].at(0);
      char chain_id1 = result_parts[1].at(0);
      char resn1 = result_parts[2].at(0);
      int resi1 = atoi(result_parts[3].c_str());
      char chain_id2 = result_parts[4].at(0);
      char resn2 = result_parts[5].at(0);
      int resi2 = atoi(result_parts[6].c_str());
      int cluster_id = atoi(result_parts[7].c_str());

      // preberemo rezultate samo, ce so za klaster z ugodnimi score-i
      if (cluster.find(cluster_id) != cluster.end()) {
	char ch1 = chain_id1;
	int rei1 = resi1;
	char ch2 = chain_id2;
	int rei2 = resi2;
	char ren2 = resn2;
	// ce sta proteina zamenjana, potem je ravno obratno
	if (swap) {
	  ch1 = chain_id2;
	  rei1 = resi2;
	  ch2 = chain_id1;
	  rei2 = resi1;
	  ren2 = resn1;
	}

	// pointer na alignan protein Molecule* dolocimo kasneje
	aligned.insert(make_pair(cluster_id, make_pair(new Residue(ch1, rei1), new Residue( ch2, ren2, rei2, cluster[cluster_id], NULL, (flx == 'Y') ? true : false ))));
	
	/* ker uporabljamo multimap, v katerem se klastri lahko prekrivajo, moramo dolociti dolzino alignmenta posebej */
	rid.insert(ChResi(ch1, rei1));
    }
  }
  unique_seq_size = rid.size();
  return true;
}

