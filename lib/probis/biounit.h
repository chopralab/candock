#ifndef _BIOUNIT_H
#define _BIOUNIT_H


class BioUnit {
 public:
  ~BioUnit() {
    /* sprostimo pomnilnik za rotacijske matrike in translacijske vektorje biounitov */
    for (map<int, gsl_matrix*>::iterator it = U.begin(); it != U.end(); it++) {
      gsl_matrix_free(it->second);
    }
    for (map<int, gsl_vector*>::iterator it = t.begin(); it != t.end(); it++)
      gsl_vector_free(it->second);
  }

  map<int, gsl_matrix*> U;
  map<int, gsl_vector*> t;
//  set<char> chain_id;
  multimap<int, char> chain_id;
  vector<string> rota_pdb; // pdb file v tekstovni obliki
  string filename;

  bool najdi(string c_id) {
//    for (set<char>::iterator it = chain_id.begin(); it != chain_id.end(); it++) 
    for (multimap<int,char>::iterator it = chain_id.begin(); it != chain_id.end(); it++) 
      if (c_id.find(it->second) != string::npos) 
        return true; 
    return false; 
  }

  bool is_in_rota_num(int rota_num, char c_id) {
    pair<multimap<int,char>::iterator,multimap<int,char>::iterator> ret = chain_id.equal_range(rota_num);
    for (multimap<int,char>::iterator it = ret.first; it != ret.second; it++) 
      if (c_id == it->second) 
        return true; 
    return false; 
  }

};




#endif // _BIOUNIT_H
