/*  Multi-neighborhood tabu search for the maximum weight clique problem
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * This program demonstrates the use of the heuristic algorithm for solving
 * the problem of the maximum weight clique problem. For more information about this algorithm, 
 * please visit the website http://www.info.univ-angers.fr/pub/hao/mnts_clique.html or email to: wu@info-univ.angers.fr.
 */             
#ifndef MNTS_H
#define MNTS_H
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <exception>
#include <vector>
#include <memory>
#include <string.h>
#include "helper/error.hpp"
#include "helper/debug.hpp"
#include "helper/inout.hpp"
//~ #include "graph2.hpp"

using namespace std;

namespace Glib {
	typedef vector<vector<bool>> AdjacencyMatrix;
};

class MNTS {
	Glib::AdjacencyMatrix &Edge;   // adjacent matrix
	unique_ptr<int[]> vectex;
	unique_ptr<int[]> funch;
	unique_ptr<int[]> address;
	unique_ptr<int[]> tabuin;
	int Max_Vtx,  Max_Iter; 
	int f;
	int fbest;
	//~ const unique_ptr<int[]> &adjaclen;
	//~ const unique_ptr<unique_ptr<int[]>[]> &adjacMatrix;
	unique_ptr<int[]> adjaclen;
	unique_ptr<unique_ptr<int[]>[]> adjacMatrix;
	unique_ptr<int[]> cruset;
	int len;
	//~ int tm1;
	int tm2;
	unique_ptr<int[]> C0;
	unique_ptr<int[]> C1;
	const int *We;
	//~ const vector<int> &We;
	unique_ptr<int[]> BC;
	int len0;
	int len1;
	unique_ptr<int[]> TC1;
	int Iter;
	int TABUL;
	int Wf;
	int Wbest;
	unique_ptr<int[]> FC1;
	//~ unique_ptr<int[]> Tbest;
	vector<vector<int>> &qmax;
	//~ unique_ptr<int[]> TTbest;
	int Waim;
	int Titer;
	int len_best;
	int len_W;
	unique_ptr<int[]> Iteration;
	unique_ptr<int[]> len_used;
	unique_ptr<int[]> W_used;
	int len_improve;
	int len_time;
	//~ int Wmode;
	//int TABUL0 = 5;
	int iter; // number of iterations
	//~ void initialize(bool **);
	void initialize();
	int randomInt( int n ) { return rand() % n; }
	void clearGamma();
	int selectC0();
	int WselectC0();
	int expand(int);
	int selectC1();
	int WselectC1();
	int plateau( int );
	int Mumi_Weigt();
	int backtract();
	int tabu( int );
	//~ void verify();
	void output();
	void max_tabu(int);
public:
	MNTS(vector<vector<int>> &qm, Glib::AdjacencyMatrix &conn, const int *weight, const int ii=300, const int w=100, 
		const int lni=10) : qmax(qm), Edge(conn), Max_Vtx(conn.size()), We(weight), 
		iter(ii), Waim(w), len_improve(lni), len_time(int (100000000 / lni) + 1), TABUL(7), len_best(0) {
#ifndef NDEBUG
		stringstream ss;
		for (int i=0; i < Max_Vtx; i++) ss << We[i] << " ";
		dbgmsg("weight = " << ss.str());
#endif
		initialize();
		dbgmsg("after initialize num. of iter = " << iter);
		for(int i = 0; i < iter; i++ ) {
			dbgmsg("doing iter " << i);
			max_tabu(i);
		}
#ifndef NDEBUG
		output();
#endif
	}
	//~ MNTS(vector<vector<int>> &qm, Glib::AdjacencyMatrix &conn, const int *weight, 
		//~ const unique_ptr<unique_ptr<int[]>[]> &adjacMatrix, 
		//~ const unique_ptr<int[]> &adjaclen, const int ii=300, const int w=600, 
		//~ const int lni=10) : qmax(qm), Edge(conn), Max_Vtx(conn.size()), We(weight), 
		//~ adjacMatrix(adjacMatrix), adjaclen(adjaclen), iter(ii), Waim(w), 
		//~ len_improve(lni), len_time(int (100000000 / lni) + 1), TABUL(7), len_best(0) {
//~ #ifndef NDEBUG
		//~ stringstream ss;
		//~ for (int i=0; i < Max_Vtx; i++) ss << We[i] << " ";
		//~ dbgmsg("weight = " << ss.str());
		//~ for (int i=0; i < Max_Vtx; i++) {
			//~ for (int j=0; j < Max_Vtx; j++) {
				//~ dbgmsg("conn[" << i << "][" << j << "] = " 
					//~ << boolalpha << conn[i][j]);
			//~ }
		//~ }
		//~ for (int i=0; i < Max_Vtx; i++) {
			//~ for (int j=0; j < adjaclen[i]; j++) {
				//~ dbgmsg("adjaclen[" << i << "] = " << adjaclen[i] << " adjacMatrix["
					//~ << i << "][" << j << "] = " << adjacMatrix[i][j]);
			//~ }
		//~ }
//~ #endif
		//~ initialize();
		//~ dbgmsg("after initialize num. of iter = " << iter);
		//~ for(int i = 0; i < iter; i++ ) {
			//~ dbgmsg("doing iter " << i);
			//~ max_tabu(i);
		//~ }
//~ #ifndef NDEBUG
		//~ output();
//~ #endif
	//~ }
};
#endif
