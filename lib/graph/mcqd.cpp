/*
Copyright (c) 2007-2019, Janez Konc
All rights reserved.
If you use this program, please cite: 
Janez Konc and Dusanka Janezic. An improved branch and bound algorithm for the 
maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.
More information at: http://insilab.org
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "mcqd.hpp"
#include "helper/help.hpp"
#include <utility>

//~ Maxclique::Maxclique (const Array2d<bool> &conn, const vector<double> &scores, const float tt) : __conn(conn), __scores(scores), pk(0), level(1), Tlimit(tt), V(conn.get_szi()), Q(conn.get_szi()), QMAX(conn.get_szi()) {
Maxclique::Maxclique (const Array2d<bool> &conn, const float tt) : __conn(conn), pk(0), level(1), Tlimit(tt), V(conn.get_szi()), QMAX(conn.get_szi()), Q(conn.get_szi()) {
  if (conn.get_szi()==0) throw Error("WARNING: Graph is empty\n(ussually means that the molecule is built of non-seed fragments only and you need to increase the binding site size)."); // fixes issue #116
  for (size_t i=0; i < conn.get_szi(); i++) V.push(i);
  C = new ColorClass[conn.get_szi() + 1];
  for (size_t i=0; i < conn.get_szi() + 1; i++) C[i].init(conn.get_szi() + 1);
  S = new StepCount[conn.get_szi() + 1];
}

vector<vector<unsigned short int>>& Maxclique::mcq(const int minsz){ 

  // search only for cliques of size >= minsz
  for (int i = 0; i < minsz - 1; ++i) QMAX.push(0);

  V.set_degrees(*this);
  V.sort();
  V.init_colors();
  expand(V);

  //help::memusage("after expand");

  return QMAXES;
}

vector<vector<unsigned short int>>& Maxclique::mcqdyn(const int minsz){ 

  // search only for cliques of size >= minsz
  for (int i = 0; i < minsz - 1; ++i) QMAX.push(0);

  V.set_degrees(*this);
  V.sort();
  V.init_colors();
  for (int i=0; i < V.size() + 1; i++) {
    S[i].set_i1(0);
    S[i].set_i2(0);
  }
  expand_dyn(V);

  //help::memusage("after expand");

  return QMAXES;
}

//~ set<pair<vector<unsigned short int>, double>, Maxclique::comp>&  Maxclique::mcqw(const int minsz) { 
//~ 
  //~ // search only for cliques of size >= minsz
  //~ for (int i = 0; i < minsz - 1; ++i) QMAX.push(0);
//~ 
  //~ V.set_degrees(*this);
  //~ V.sort();
  //~ V.init_colors();
  //~ expand_weight(V);
//~ 
  //~ //help::memusage("after expand");
//~ 
  //~ return QMAXES_WEIGHT;
//~ }

void Maxclique::Vertices::init_colors() { 
  const int max_degree = v[0].get_degree();
  for (int i = 0; i < max_degree; i++)
    v[i].set_degree(i + 1);
  for (int i = max_degree; i < sz; i++)
    v[i].set_degree(max_degree + 1);
}

void Maxclique::Vertices::set_degrees(Maxclique &m) { 
  for (int i=0; i < sz; i++) {
    int d = 0;
    for (int j=0; j < sz; j++)
      if (m.connection(v[i].get_i(), v[j].get_i())) d++;
    v[i].set_degree(d);
  }
}

bool Maxclique::cut1(const int pi, const ColorClass &A) {
  for (int i = 0; i < A.size(); i++)
    if (connection(pi, A.at(i)))
      return true;
  return false;
}

void Maxclique::cut2(const Vertices &A, Vertices &B) {
  for (int i = 0; i < A.size() - 1; i++) {
    if (connection(A.end().get_i(), A.at(i).get_i()))
      B.push(A.at(i).get_i());
  }
}

void Maxclique::color_sort(Vertices &R) {
  int j = 0;
  int maxno = 1;
  int min_k = QMAX.size() - Q.size() + 1;
  C[1].rewind();
  C[2].rewind();
  int k = 1;
  for (int i=0; i < R.size(); i++) {
    int pi = R.at(i).get_i();
    k = 1;
    while (cut1(pi, C[k]))
      k++;
    if (k > maxno) {
      maxno = k;
      C[maxno + 1].rewind();
    }
    C[k].push(pi);
    if (k < min_k) {
      R.at(j++).set_i(pi);
    }
  }
  if (j > 0) R.at(j-1).set_degree(0);
  if (min_k <= 0) min_k = 1;
  for (k = min_k; k <= maxno; k++)
    for (int i = 0; i < C[k].size(); i++) {
      R.at(j).set_i(C[k].at(i));
      R.at(j++).set_degree(k);
    }
}

void Maxclique::expand(Vertices R) {
  if (QMAXES.size() > 10000000) return;
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      //~ cout << "QMs = " << QMAX.size() << endl;
      //~ cout << "Qs = " << Q.size() << endl;
      //~ if (Rp.size()) {
      if (Rp.size() && Q.size() < QMAX.size() + 1) {
        color_sort(Rp);
		pk++;
        expand(Rp);
      }
      else if (Q.size() > QMAX.size()) { 
        //~ cout << "step = " << pk << " Q.size() = " << Q.size() << endl; 
		//~ QMAXES.push_back(Q); //akj

		vector<unsigned short int> lastQ;
        for (int i=0; i < Q.size(); ++i) { 
			lastQ.push_back(static_cast<unsigned short int>(Q.at(i))); 
		}
        QMAXES.push_back(lastQ); //akj
      }    
      Rp.dispose();
      Q.pop();
    }
    else {
      return;
    }
    R.pop();
  }
}

//~ void Maxclique::expand_weight(Vertices R) {
  //~ if (QMAXES_WEIGHT.size() > 100000) return;
  //~ while (R.size()) {
    //~ if (Q.size() + R.end().get_degree() > QMAX.size()) {
      //~ Q.push(R.end().get_i());
      //~ Vertices Rp(R.size());
      //~ cut2(R, Rp);
      //~ if (Rp.size()) {
        //~ color_sort(Rp);
		//~ pk++;
        //~ expand(Rp);
      //~ }
      //~ else if (Q.size() > QMAX.size()) { 
		//~ vector<unsigned short int> lastQ;
		//~ double energy = 0.0;
        //~ for (int i=0; i < Q.size(); ++i) { 
			//~ lastQ.push_back(static_cast<unsigned short int>(Q.at(i)));
			//~ energy += __scores[Q.at(i)];
		//~ }
		//~ QMAXES_WEIGHT.insert(make_pair(lastQ, energy));
		//~ if (QMAXES_WEIGHT.size() > 100)
			//~ QMAXES_WEIGHT.erase(--QMAXES_WEIGHT.end());
      //~ }    
      //~ Rp.dispose();
      //~ Q.pop();
    //~ }
    //~ else {
      //~ return;
    //~ }
    //~ R.pop();
  //~ }
//~ }

void Maxclique::expand_dyn(Vertices R) {
  S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
  S[level].set_i2(S[level - 1].get_i1());
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      //~ if (Rp.size()) {
      if (Rp.size() && Q.size() < QMAX.size() + 1) {
        if ((float)S[level].get_i1()/++pk < Tlimit) {
          degree_sort(Rp);
        }
		color_sort(Rp);
		S[level].inc_i1();
		level++;
		expand_dyn(Rp);
		level--;
      }
      else if (Q.size() > QMAX.size()) { 
		//~ cout << "step = " << pk << " current max. clique size = " << Q.size() << endl; 
		//~ QMAXES.push_back(Q);
		vector<unsigned short int> lastQ;
        for (int i=0; i < Q.size(); ++i) { 
			lastQ.push_back(static_cast<unsigned short int>(Q.at(i))); 
		}
        QMAXES.push_back(lastQ); //akj
      }    
      Rp.dispose();
      Q.pop();
    }
    else {
      return;
    }
    R.pop();
  }
}

