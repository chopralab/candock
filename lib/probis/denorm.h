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

#ifndef _DENORM_H
#define _DENORM_H

#include "const.h"

class ClusterData;
class Residue;

class Denorm {
  map<int, ClusterData*> cluster;
  map<int, ClusterData*>::iterator curr_cl_it;
  multimap<int, pair<Residue*, Residue*> > aligned;
  multimap<int, pair<Residue*, Residue*> >::iterator curr_al_it;
  pair<multimap<int, pair<Residue*, Residue*> >::iterator, multimap<int, pair<Residue*, Residue*> >::iterator > equal_al;
  string pdb_id;
  string chain_id;
  int unique_seq_size;
 public:
  void free(int);
  bool read(string&, bool=true);
  bool read_one(string&);
  gsl_matrix* get_U(int);
  gsl_vector* get_t(int);
  vector<pair<ChResi, int> > get_cons(int);
  string get_pdb_id() { return pdb_id; }
  string get_chain_id() { return chain_id; }
  ClusterData* get_next_cluster() { if (curr_cl_it == cluster.end()) { return NULL; } else { ClusterData *res = curr_cl_it->second; curr_cl_it++; return res; } }
  void reset_cluster() { curr_cl_it = cluster.begin(); }
  bool get_next_residue(Residue *&r1, Residue *&r2) { if (curr_al_it != equal_al.second) { r1 = curr_al_it->second.first; r2 = curr_al_it->second.second; curr_al_it++; return true;} else return false; }
  void reset_residue(int cluster_id) { equal_al = aligned.equal_range(cluster_id); curr_al_it = equal_al.first; }
  int get_unique_seq_size() { return unique_seq_size; }
};

#endif
