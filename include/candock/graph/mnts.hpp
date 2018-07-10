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
#include "candock/helper/error.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/inout.hpp"
//~ #include "graph2.hpp"

namespace candock {
namespace graph {

typedef std::vector<std::vector<bool>> AdjacencyMatrix;

class MNTS {
	AdjacencyMatrix &Edge;   // adjacent matrix
	std::unique_ptr<int[]> vectex;
	std::unique_ptr<int[]> funch;
	std::unique_ptr<int[]> address;
	std::unique_ptr<int[]> tabuin;
	int Max_Vtx,  Max_Iter; 
	int f;
	int fbest;
	//~ const unique_ptr<int[]> &adjaclen;
	//~ const unique_ptr<unique_ptr<int[]>[]> &adjacMatrix;
	std::unique_ptr<int[]> adjaclen;
	std::unique_ptr<std::unique_ptr<int[]>[]> adjacMatrix;
	std::unique_ptr<int[]> cruset;
	int len;
	//~ int tm1;
	int tm2;
	std::unique_ptr<int[]> C0;
	std::unique_ptr<int[]> C1;
	const int *We;
	//~ const vector<int> &We;
	std::unique_ptr<int[]> BC;
	int len0;
	int len1;
	std::unique_ptr<int[]> TC1;
	int Iter;
	int TABUL;
	int Wf;
	int Wbest;
	std::unique_ptr<int[]> FC1;
	//~ unique_ptr<int[]> Tbest;
	std::vector<std::vector<int>> &qmax;
	//~ unique_ptr<int[]> TTbest;
	int Waim;
	int Titer;
	int len_best;
	int len_W;
	std::unique_ptr<int[]> Iteration;
	std::unique_ptr<int[]> len_used;
	std::unique_ptr<int[]> W_used;
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
	MNTS(std::vector<std::vector<int>> &qm, AdjacencyMatrix &conn, const int *weight, const int ii=300, const int w=100, 
		const int lni=10) : Edge(conn), Max_Vtx(conn.size()), We(weight), TABUL(7),
                qmax(qm), Waim(w), len_best(0), len_improve(lni), len_time(int (100000000 / lni) + 1), iter(ii) {
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

};
}
}

#endif
