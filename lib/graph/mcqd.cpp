#include "mcqd.hpp"
#include "helper/help.hpp"

Maxclique::Maxclique (Array2d<bool> &conn, const float tt) : e(conn), pk(0), level(1), Tlimit(tt), V(conn.get_szi()), Q(conn.get_szi()), QMAX(conn.get_szi()) {
  if (conn.get_szi()==0) throw exc_empty;
  for (int i=0; i < conn.get_szi(); i++) V.push(i);
  C = new ColorClass[conn.get_szi() + 1];
  for (int i=0; i < conn.get_szi() + 1; i++) C[i].init(conn.get_szi() + 1);
  S = new StepCount[conn.get_szi() + 1];
}

std::vector<std::vector<short int>>& Maxclique::_mcq(const int minsz, bool dyn){ 

  // search only for cliques of size >= minsz
  for (int i = 0; i < minsz - 1; ++i) QMAX.push(0);

  V.set_degrees(*this);
  V.sort();
  V.init_colors();
  if (dyn) {
    for (int i=0; i < V.size() + 1; i++) {
      S[i].set_i1(0);
      S[i].set_i2(0);
    }
    expand_dyn(V);
  }
  else
    expand(V);

  help::memusage("after expand");

  QMAXES.shrink_to_fit();
  
  help::memusage("after shrink to fit");

  return QMAXES;
}

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
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      //~ std::cout << "QMs = " << QMAX.size() << std::endl;
      //~ std::cout << "Qs = " << Q.size() << std::endl;
      if (Rp.size()) {
        color_sort(Rp);
		pk++;
        expand(Rp);
      }
      else if (Q.size() > QMAX.size()) { 
        //~ std::cout << "step = " << pk << " Q.size() = " << Q.size() << std::endl; 
		//~ QMAXES.push_back(Q); //akj
		std::vector<short int> lastQ;
        for (int i=0; i < Q.size(); ++i) { 
			lastQ.push_back(static_cast<short int>(Q.at(i))); 
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

void Maxclique::expand_dyn(Vertices R) {
  S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
  S[level].set_i2(S[level - 1].get_i1());
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
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
		//~ std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
		//~ QMAXES.push_back(Q);
		std::vector<short int> lastQ;
        for (int i=0; i < Q.size(); ++i) { 
			lastQ.push_back(static_cast<short int>(Q.at(i))); 
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

