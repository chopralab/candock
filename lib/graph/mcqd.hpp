/*
    Copyright 2007-2012 Janez Konc 

    If you use this program, please cite: 
    Janez Konc and Dusanka Janezic. An improved branch and bound algorithm for the 
    maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.

    More information at: http://www.sicmm.org/~konc

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MCQD_H
#define MCQD_H

#include <iostream>
#include <algorithm>
#include <assert.h>
#include <exception>
#include <vector>

#ifdef DBG
using namespace std;
#endif

class Maxclique {
  const bool* const* e;
  int pk, level;
  const float Tlimit;
  class Vertices {
    class Vertex {
      int i, d;
    public:
      void set_i(const int ii)  { i = ii; }
      int get_i() const { return i; }
      void set_degree(int dd) { d = dd; }
      int get_degree() const { return d; }
    };
    Vertex *v;
    int sz;
    static bool desc_degree(const Vertex vi, const Vertex vj) { return (vi.get_degree() > vj.get_degree()); }
  public:
#ifdef DBG
    void dbg_v(const string msg="") const {
    //  std::cout << msg << " Vertices: [";
      for (int i=0; i < sz; i++) 
//	std::cout << "(" << v[i].get_i() << "," << v[i].get_degree() << ") ";
  //    std::cout << "]" << std::endl;
    }
#endif
    Vertices(int size) : sz(0) { v = new Vertex[size]; }
    ~Vertices () {}
    void dispose() { if (v) delete [] v; }
    void sort() { std::sort(v, v+sz, desc_degree); }
    void init_colors();
    void set_degrees(Maxclique&);
    int size() const { return sz; }
    void push(const int ii) { v[sz++].set_i(ii); };
    void pop() { sz--; };
    Vertex& at(const int ii) const { return v[ii]; };
    Vertex& end() const { return v[sz - 1]; };
  };
  class ColorClass {
    int *i;
    int sz;
  public:
#ifdef DBG
    void dbg_i(const string msg="") const {
      //std::cout << msg << " Class: [";
      for (int ii=0; ii < sz; ii++) 
	//std::cout << i[ii] << " ";
      //std::cout << "]" << std::endl;
    }
#endif
    ColorClass() : sz(0), i(0) {}
    ColorClass(const int sz) : sz(sz), i(0) { init(sz); }
    ~ColorClass() { if (i) delete [] i;
    }
    ColorClass(const ColorClass& dh) { // copy
		init(dh.sz);
      for (int j = 0; j < dh.sz; j++) i[j] = dh.i[j];
      sz = dh.sz;
    }
    void init(const int sz) { i = new int[sz]; rewind(); }
    void push(const int ii) { i[sz++] = ii; };
    void pop() { sz--; };
    void rewind() { sz = 0; };
    int size() const { return sz; }
    int& at(const int ii) const { return i[ii]; }
    ColorClass& operator=(const ColorClass& dh) {
      for (int j = 0; j < dh.sz; j++) i[j] = dh.i[j];
      sz = dh.sz;
      return *this;
    }
  };
  Vertices V;
  ColorClass *C, QMAX, Q;
  std::vector<ColorClass> QMAXES;
  class StepCount {
    int i1, i2;
  public:
    StepCount() : i1(0), i2(0) {}
    void set_i1(const int ii)  { i1 = ii; }
    int get_i1() const { return i1; }
    void set_i2(const int ii)  { i2 = ii; }
    int get_i2() const { return i2; }
    void inc_i1()  { i1++; }
  };
  StepCount *S;
  bool connection(const int i, const int j) const { return e[i][j]; }
  bool cut1(const int, const ColorClass&);
  void cut2(const Vertices&, Vertices&);
  void color_sort(Vertices&);
  void expand(Vertices);
  void expand_dyn(Vertices);
  void _mcq(std::vector<std::vector<int>>&, int&, bool);
  //~ void _mcq(int**&, int&, bool);
  //~ void _mcq(int*&, int&, bool);
  void degree_sort(Vertices &R) { R.set_degrees(*this); R.sort(); }
public:
#ifdef DBG
  void dbg_C() const {
    for (int i=0; i < V.size(); i++) {
      //std::cout << "C["<< i << "] : ";
      C[i].dbg_i();
    }
  }
  void dbg_conn() const {
    for (int i=0; i < V.size(); i++) {
      for (int j=0; j < V.size(); j++) {
	//std::cout <<e[i][j];
      }
    //  std::cout<< std::endl;
    }
  }
#endif
  class myexception : public std::exception {
    const char* what() const throw() { return "WARNING: Graph is empty."; }
  } exc_empty;
  Maxclique(const bool* const*, const int, const float=0.025);
  int steps() const { return pk; }
  void mcq(std::vector<std::vector<int>> &maxclique, int &sz) { _mcq(maxclique, sz, false); }
  //~ void mcq(int** &maxclique, int &sz) { _mcq(maxclique, sz, false); }
  void mcqdyn(std::vector<std::vector<int>> &maxclique, int &sz) { _mcq(maxclique, sz, true); }
  //~ void mcqdyn(int** &maxclique, int &sz) { _mcq(maxclique, sz, true); }
  //~ void mcqdyn(int* &maxclique, int &sz) { _mcq(maxclique, sz, true); }
  ~Maxclique() {
    //std::cout << "Maxclique::~Maxclique" << std::endl;
    if (C) delete [] C;
    if (S) delete [] S;
    V.dispose();
  };
};
#endif
