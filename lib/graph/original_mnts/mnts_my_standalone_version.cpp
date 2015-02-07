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
//~ #infdef MNTS_H
//~ #define MNTS_H
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <exception>
#include <vector>
#include <string.h>
#include "helper/error.hpp"
#include "helper/debug.hpp"
#include "helper/inout.hpp"
//~ #include "pdbreader/benchmark.hpp"

using namespace std;

namespace Glib {
	typedef vector<vector<bool>> AdjacencyMatrix;
};

class MNTS {
	string file_name;
	Glib::AdjacencyMatrix &Edge;   // adjacent matrix
	//~ int **Edge;   // adjacent matrix
	int *vectex;
	int *funch;
	int *address;
	int *tabuin;
	int Max_Vtx,  Max_Iter; 
	int f;
	int fbest;
	//~ int *adjaclen;
	//~ int **adjacMatrix;
	unique_ptr<int[]> &adjaclen;
	unique_ptr<unique_ptr<int[]>[]> &adjacMatrix;
	int *cruset;
	int len;
	//~ int tm1;
	int tm2;
	int *C0;
	int *C1;
	int *We;
	int *BC;
	int len0;
	int len1;
	int *TC1;
	int Iter;
	int TABUL;
	int Wf;
	int Wbest;
	int *FC1;
	int *Tbest;
	vector<vector<int>> qmax;
	int *TTbest;
	int Waim;
	int Titer;
	//~ double starting_time, finishing_time, avg_time;
	int len_best;
	//~ int len_best = 0;
	int len_W;
	int Iteration[ 100 ];
	//~ double time_used[ 100 ];
	int len_used[ 100 ];
	int W_used[ 100 ];
	//~ char outfilename[30];
	int len_improve;
	int len_time;
	//~ int Wmode;
	//int TABUL0 = 5;
public:
	MNTS(const string &fn, const int w, const int lni) : file_name(fn), Waim(w), len_improve(lni), len_time(int (100000000 / lni) + 1), TABUL(7), len_best(0) {
		initialize();
	}
	//~ MNTS(std::vector<std::vector<int>> &maxclique, , const int w, const int lni
	// section 0, initiaze
	void initialize() {
		//~ Benchmark::reset();
		srand( (unsigned) time( NULL ) );
		vector<string> dimacs;
		Inout::read_file(file_name, dimacs);
		int max_edg=0;
		int nb_edg;
		for (auto &line : dimacs) {
			if (line.at(0) == 'p') {
				stringstream ss;
				ss << line;
				char ctmp;
				ss >> ctmp >> Max_Vtx >> nb_edg;
				dbgmsg("Number of vertexes = " << Max_Vtx);
				dbgmsg("Number of edges = " << nb_edg);
				vectex = new int [Max_Vtx];
				funch  = new int [Max_Vtx];
				address= new int [Max_Vtx];
				tabuin = new int [Max_Vtx];
				//~ adjaclen= new int[Max_Vtx];
				C0 = new int[Max_Vtx];
				C1 = new int[Max_Vtx];
				TC1= new int[Max_Vtx];
				We = new int[Max_Vtx];
				BC = new int[Max_Vtx];
				FC1= new int[Max_Vtx];
				//~ adjacMatrix = new int*[Max_Vtx];
				cruset = new int [2000];
				Tbest = new int[Max_Vtx];
				TTbest = new int[Max_Vtx];
				for( int x = 0 ; x < Max_Vtx ; x++ ) {
				   vectex[x]  = x;
				   address[x] = x;
				}

				Edge.resize(Max_Vtx); 
				for (auto &row : Edge) 
					row.resize(Max_Vtx, false); 
				//~ Edge = new int*[ Max_Vtx ];  
				//~ for (int x = 0 ; x < Max_Vtx ; x++ ) {
				  //~ Edge[x] = new int[Max_Vtx];
				  //~ adjacMatrix[x] = new int[Max_Vtx];
				//~ }
				for ( int x=0; x<Max_Vtx; x++ )
					for ( int y=0; y<Max_Vtx; y++ )
						Edge[x][y] = 0;
				
				
				unique_ptr<int[]> adjaclen(new int[Max_Vtx]);
				unique_ptr<unique_ptr<int[]>[]> adjacMatrix(new unique_ptr<int[]>[Max_Vtx]);
				dbgmsg("max_weight_clique --- start filling : Size of graph = " 
					<< Max_Vtx);
				for (auto &v : *this) {
					const int &x = v.get_index();
					adjaclen[x] = 0;
					adjacMatrix[x] = unique_ptr<int[]>(new int[v.size()]);
					dbgmsg("v.size() = " << v.size());
					for (auto &adj_v : v) {
						const int &y = adj_v.get_index();
						adjacMatrix[x][adjaclen[x]] = y;
						dbgmsg("max_weight_clique : x = " << x << " y = " << y
							<< " adjacMatrix = " << adjacMatrix[x][adjaclen[x]]);
						adjaclen[x]++;
					}
				}


			}
			if (line.at(0) == 'e') {
				stringstream ss;
				ss << line;
				char ctmp;
				int x1, x2;
				ss >> ctmp >> x1 >> x2;
				x1--; x2--;
				if ( x1<0 || x2<0 || x1>=Max_Vtx || x2 >=Max_Vtx ) throw Error("### Error of node : x1=" + to_string(x1) + ", x2=" + to_string(x2));
				 Edge[x1][x2]=Edge[x2][x1]=1;
				 max_edg++;	                 
			}
			if (line.at(0) == 'w') {
				stringstream ss;
				ss << line;
				char ctmp;
				int x, w;
				ss >> ctmp >> x >> w;
				x--;
				We[ x ] = w;
				dbgmsg("weight["<<x<<"] = " << w);
				//~ We[ x ] = (x+1)%Wmode + 1;
				BC[ x ] = 0;
			}
		}
		dbgmsg("Density = " << (float) max_edg/(Max_Vtx*(Max_Vtx-1)));
		if ( 0 && max_edg != nb_edg ) throw Error("### Error de lecture du graphe, nbre aretes : annonce=" + to_string(nb_edg) + ", lu=" + to_string(max_edg));
		for( int x=0 ; x<Max_Vtx; x++ )
			for( int y=0 ; y<Max_Vtx; y++ )
				Edge[x][y] = 1 - Edge[x][y];
		for( int x=0 ; x<Max_Vtx; x++ )
			Edge[x][x] = 0;
		for( int x=0 ; x<Max_Vtx; x++ ) {
			adjaclen[x] = 0;
			for( int y=0; y<Max_Vtx; y++ ) {
				if( Edge[x][y] == 1 ) {
					adjacMatrix[x][adjaclen[x]] = y;
					adjaclen[x]++;
				}
			}
		}
		//~ cout << "time to read data " << Benchmark::seconds_from_start() << " wallclock seconds" << endl;
	}
	int randomInt( int n ) {
		return rand() % n;
	}
	void clearGamma() {
		int i, j, k, l;
		int tm1 = Max_Vtx*sizeof( int );
		memset( vectex, 0, tm1 );
		memset(  funch, 0, tm1 );
		memset(address, 0, tm1 );
		memset( tabuin, 0, tm1 );
		for( i = 0; i < Max_Vtx; i++ ) {
			C0[ i ] = i;
			address[ i ] = i;
		}
		len0 = Max_Vtx;
		len1 = 0;
		len = 0;
		Wf = 0;
		Wbest = 0;
	}
	int selectC0() {
		int i, j, k, l, m;
		l = 0;
		if( len0 > 30 ) {
			k = randomInt( len0 );
			return k;
		}
		for( i = 0; i < len0; i++ ) {
			k = C0[ i ];
			if( tabuin[ k ] <= Iter )
				TC1[ l++ ] = i;
		}
		if( l == 0 )
		return -1;
		else {
			k = randomInt( l );
			k = TC1[ k ];
			return k;
		}
	}
	int WselectC0( ) {
		int i, j, k, l1, l2, w1, w2, m;
		l1 = 0;
		l2 = 0;
		w1 = 0;
		w2 = 0;
		for( i = 0; i < len0; i++ ) {
			k = C0[ i ];
			if( tabuin[ k ] <= Iter ) {
				if( We[ k ] > w1 ) {
					l1 = 0;
					w1 = We[ k ];
					FC1[ l1++ ] = i;
				}
				else if ( We[ k ] >= w1 ) {
					FC1[ l1++ ] = i;
				}
			} else {
				if( We[ k ] > w2 ) {
					l2 = 0;
					w2 = We[ k ];
					TC1[ l2++ ] = i;
				}
				else if ( We[ k ] >= w2 ) {
					TC1[ l2++ ] = i;
				}
			}
		}
		if( (l2 > 0) && ( w2 > w1 ) && ((w2+Wf)>Wbest) ) {
			k = randomInt( l2 );
			k = TC1[ k ];
			//cout << "yes in aspiration w2+Wf = " << w2+Wf << endl;
			//getchar();
			return k;
		}  
		else if( l1 > 0 ) {
			k = randomInt( l1 );
			k = FC1[ k ];
			return k;
		} else {
			return -1;
		}
	}
	int expand(int SelN)
	{
	    int i, j, k, k1, l, am, m, n, n1;
	    
	    m = C0[ SelN ];
	    cruset[ len++ ] = m;
	    vectex[ m ] = 1;
	    Wf = Wf + We[ m ];
	    
	    len0--;
	    n1 = C0[ len0 ];
	    k1 = address[ m ];
	    C0[ k1 ] = n1;
	    address[ n1 ] = k1;
	    
	    for( i = 0; i < adjaclen[ m ]; i++ )
	    {
	       n = adjacMatrix[ m ][ i ];
	       funch[ n ]++;
	       if( funch[ n ] == 1 )
	       {
	           k1 = address[ n ];
	           len0--;
	           n1 = C0[ len0 ];
	           C0[ k1 ] = n1;
	           address[ n1 ] = k1;
	           
	           C1[ len1 ] = n;
	           address[ n ] = len1;
	           len1++;
	           BC[ n ] = m;
	       }
	       else if( funch[ n ] == 2 )
	       {
	           len1--;
	           n1 = C1[ len1 ];
	           k1 = address[ n ];
	           C1[ k1 ] = n1;
	           address[ n1 ] = k1;
	       }
	    } 
	    
	    if( Wf > Wbest )
	     {
	        Wbest = Wf;
	        len_best = len;
			//~ dbgmsg("max.weight clique vertices = ");
			//~ qmax.back().clear();
			//~ for (i = 0; i < len; i++) {
				//~ qmax.back().push_back(cruset[i]);
			//~ }
	        //~ for( i = 0; i < Max_Vtx; i++ ) {
				//~ if (vectex[i] > 0) {
				//~ if (vectex[i] == 1) {
					//~ qmax.back().push_back(i);
					//~ dbgmsg("vectex[" << i << "] = " << vectex[i]);
				//~ }
	            //~ Tbest[ i ] = vectex[ i ];
	        //~ }
	        //~ verify();
	        /*for( i = 0; i < Max_Vtx; i++ )
	        {
	            Tbest[ i ] = vectex[ i ];
	        }*/
	     }
	    
	    return 1;   
	}
	
	int selectC1( )
	{
	    int i, j, k, l, m;
	    l = 0;
	    for( i = 0; i < len1; i++ )
	    {
	       k = C1[ i ];
	       if( tabuin[ k ] <= Iter )
	         TC1[ l++ ] = i;
	    }
	    if( l == 0 )
	      return -1;
	    else
	    {
	        k = randomInt( l );
	        k = TC1[ k ];
	        return k;
	    }
	}
	
	int WselectC1( )
	{
	     int i, j, k, l, l1, l2, wmn, w1, w2, m, n;
	     l1 = 0;
	     l2 = 0;
	     w1 = -1000000;
	     w2 = -1000000;
	     l = 0;
	     for( i = 0; i < len1; i++ )
	     {
	         m = C1[ i ];
	         n = BC[ m ];
	         if( (vectex[ n ] == 1) && (Edge[ m ][ n ] == 1) )
	           l++;
	         else
	         {
	             for( j = 0; j < len; j++ )
	             {
	                k = cruset[ j ];
	                if( Edge[ m ][ k ] == 1 )
	                  break;
	             }
	             BC[ m ] = k;
	         }
	     }
	     //cout << "len1 = " << len1 << " l = " << l << endl;
	     for( i = 0; i < len1; i++ )
	     {
	         m = C1[ i ];
	         n = BC[ m ];
	         wmn = We[ m ] - We[ n ];
	         if( tabuin[ m ] <= Iter )
	         {
	             if( wmn > w1 )
	             {
	                l1 = 0;
	                w1 = wmn;
	                FC1[ l1++ ] = i;
	             }
	             else if ( wmn >= w1 )
	             {
	                FC1[ l1++ ] = i;
	             }
	         }
	         else
	         {
	             if( wmn > w2 )
	             {
	                l2 = 0;
	                w2 = wmn;
	                TC1[ l2++ ] = i;
	             }
	             else if ( wmn >= w2 )
	             {
	                TC1[ l2++ ] = i;
	             }
	         }
	     }
	     
	     if( (l2 > 0) && ( w2 > w1 ) && ((w2+Wf)>Wbest) )
	     {
	         k = randomInt( l2 );
	         k = TC1[ k ];
	         return k;
	     }  
	     else if( l1 > 0 )
	     {
	         k = randomInt( l1 );
	         k = FC1[ k ];
	         return k;
	     }
	     else
	     {
	         return -1;
	     }
	}
	
	int plateau( int SelN )
	{
	     int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;
	     
	     m = C1[ SelN  ];
	     for(ti = 0; ti < len; ti++)
	     {
	         m1 = cruset[ ti ];
	         if( Edge[ m1 ][ m ] == 1 )
	            break;
	     }
	     
	     Wf = Wf + We[ m ] - We[ m1 ];
	     
	     //the expand process, put m into the current independent set
	     vectex[ m ] = 1;
	     cruset[ len++ ] = m;
	     //delete m from C1
	     k1 = address[ m ];
	     len1--;
	     n1 = C1[ len1 ];
	     C1[ k1 ] = n1;
	     address[ n1 ] = k1;
	     
	     for( i = 0; i < adjaclen[ m ]; i++ )
	     {
	        n = adjacMatrix[ m ][ i ];
	        funch[ n ]++;
	        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
	        {
	             //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
	             k1 = address[ n ];
	             len0--;
	             n1 = C0[ len0 ];
	             C0[ k1 ] = n1;
	             address[ n1 ] = k1;
	             
	             C1[ len1 ] = n;
	             address[ n ] = len1;
	             len1++;
	             BC[ n ] = m;
	           
	             //getchar();
	        }
	        if( funch[ n ] == 2 )
	        {
	            len1--;
	            n1 = C1[ len1 ];
	            k1 = address[ n ];
	            C1[ k1 ] = n1;
	            address[ n1 ] = k1;
	        }        
	     } 
	     
	     //the backtrack process, delete m1 from the current independent set
	     vectex[ m1 ] = 0;
	     //cout << "len1 = " << len1 << endl;
	     tabuin[ m1 ] = Iter + TABUL + randomInt( len1+2 );
	     len--;
	     cruset[ ti ] = cruset[ len ];
	     C1[ len1 ] = m1;
	     address[ m1 ] = len1;
	     len1++;
	     
	     for( i = 0; i < adjaclen[ m1 ]; i++ )
	     {
	        n = adjacMatrix[ m1 ][ i ];
	        funch[ n ]--;
	        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
	        {
	           k1 = address[ n ];           
	           len1--;
	           n1 = C1[ len1 ];
	           C1[ k1 ] = n1;
	           address[ n1 ] = k1;
	           
	           C0[ len0 ] = n;
	           address[ n ] = len0;
	           len0++;
	        }
	        else if( funch[ n ] == 1 )
	        {
	           C1[ len1 ] = n;
	           address[ n ] = len1;
	           len1++;
	        }
	     }
	     
	     if( Wf > Wbest )
	     {
	        Wbest = Wf;
	        len_best = len;
			//~ dbgmsg("max.weight clique vertices = ");
			//~ qmax.back().clear();
	        //~ for( i = 0; i < Max_Vtx; i++ ) {
				//~ if (vectex[i] > 0) {
					//~ qmax.back().push_back(i);
				//~ }
	            //~ Tbest[ i ] = vectex[ i ];
	        //~ }
	        //~ verify();
	        /*for( i = 0; i < Max_Vtx; i++ )
	        {
	            Tbest[ i ] = vectex[ i ];
	        }*/
	     }
	     return 1;   
	}
	
	
	int Mumi_Weigt()
	{
	    int i, j, k, l1, m;
	    int w1 = 5000000;
	    l1 = 0;
	    for( i = 0; i < len; i++ )
	    {
	       k = cruset[ i ];
	       if( We[ k ] < w1 )
	       {
	          l1 = 0;
	          w1 = We[ k ];
	          FC1[ l1++ ] = i;
	       }
	       else if ( We[ k ] <= w1 )
	       {
	          FC1[ l1++ ] = i;
	       }
	    }
	    
	    if( l1 == 0 )
	      return -1;
	    k = randomInt( l1 );
	    k = FC1[ k ];
	    return k;
	}
	
	int backtract()
	{
	     int i, j, k, l, m, m1, n, ti, k1, n1;
	     ti = Mumi_Weigt();
	     if( ti == -1 )
	      return -1;
	     m1 = cruset[ ti ];
	     Wf = Wf - We[ m1 ];
	     vectex[ m1 ] = 0;
	     tabuin[ m1 ] = Iter + TABUL;
	     len--;
	     cruset[ ti ] = cruset[ len ];
	     C0[ len0 ] = m1;
	     address[ m1 ] = len0;
	     len0++;
	     
	     for( i = 0; i < adjaclen[ m1 ]; i++ )
	     {
	        n = adjacMatrix[ m1 ][ i ];
	        funch[ n ]--;
	        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
	        {
	           k1 = address[ n ];           
	           len1--;
	           n1 = C1[ len1 ];
	           C1[ k1 ] = n1;
	           address[ n1 ] = k1;
	           
	           C0[ len0 ] = n;
	           address[ n ] = len0;
	           len0++;
	        }
	        else if( funch[ n ] == 1 )
	        {
	           C1[ len1 ] = n;
	           address[ n ] = len1;
	           len1++;
	        }
	     }
	}
	
	int tabu( int Max_Iter )
	{
	     int i, j, k, l, bestlen = 0, am, am1, ww, ww1, ww2, ti, m1;
	     Iter = 0;
	     clearGamma(); 
	     while( 1 )
	     {
	        am = selectC0();
	        if( am != -1 )
	        {
	            l = expand( am );
	            Iter++;
	            if( Wbest == Waim )
	               return Wbest;
	        }
	        else 
	            break;
	     }
	      
	     while( Iter < Max_Iter )
	     {
	        am = WselectC0();
	        am1 = WselectC1();
	        if( (am != -1) && (am1 != -1) )
	        {
	            ww = We[ C0[ am ] ];
	            ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
	        
	            if( ww > ww1 )
	            {
	                l = expand( am );
	                
	                Iter++;
	                if( Wbest == Waim )
	                   return Wbest;
	            }
	            else
	            {
	                l = plateau( am1 );
	                if( Wbest == Waim )
	                    return Wbest; 
	                Iter++;
	            }
	        }
	        else if( (am != -1) && (am1 == -1) )
	        {
	             l = expand( am );
	             if( Wbest == Waim )
	               return Wbest;
	                
	             Iter++;
	        }
	        else if( (am == -1) && (am1 != -1) )
	        {
	             ti = Mumi_Weigt();
	             m1 = cruset[ ti ];
	             ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
	             ww2 = - We[ m1 ];
	             if( ww1 > ww2 )
	             {
	                l = plateau( am1 );
	                if( Wbest == Waim )
	                    return Wbest; 
	                Iter++;
	             }
	             else
	             {
	                 k = backtract();
	                 if( k == -1 )
	                     return Wbest;
	                 Iter++;
	             }
	        }
	        else if( (am == -1) && (am1 == -1) )
	        {
	             k = backtract();
	             if( k == -1 )
	                return Wbest;
	             Iter++;
	        }
	             
	     }
	
	     return Wbest;
	}
	void verify()
	{
	     int i, j, k1, k2, l, m;
	     for( i = 0; i < Max_Vtx; i++ )
	     {
	          if( Tbest[ i ] == 1 )
	          //~ if( TTbest[ i ] == 1 )
	          {
	              for( j = i+1; j < Max_Vtx; j++ )
	              if( (Tbest[ j ] == 1) && ( Edge[ i ][ j ] == 1 ) )
	              //~ if( (TTbest[ j ] == 1) && ( Edge[ i ][ j ] == 1 ) )
	                  throw Error("hello there is something wrong");
	          }
	     }
	}
	void output() {
		int i , j, k, l, sum; 
		stringstream ss;
		for( i = 0; i < 100; i++ ) {
			ss << "sum = " << W_used[ i ] << ", iter = " << Iteration[ i ] << ", len = " << len_used[ i ] << ", qmax = "; 
			for (auto &v : qmax[i]) ss << v << " ";
			ss << endl;
		}
		ss << "\n\n the total information:" << endl;
		int wavg, iteravg, lenbb, success;
		wavg = iteravg = lenbb = success = 0;
		int best_v = 0;
		for( i = 0; i < 100; i++ ) {
			if( W_used[ i ] > best_v ) {
				best_v = W_used[ i ];  
				lenbb = len_used[ i ];
			}
	    }
		int count = 0;
		ss << "\n The best weight value for the maximum weighted problem is " << best_v << endl;
		for( i = 0; i < 100; i++ ) {
			wavg = wavg + W_used[ i ];
		}  
		double twavg = (double (wavg)) / 100 ; 
		for( i = 0; i < 100; i++ ) {
			if( W_used[ i ] == best_v ) {
				count++;
				iteravg = iteravg + Iteration[ i ];
			}
		}
		iteravg =  int ( (double (iteravg)) / count );
		ss << "avg_sum = " << twavg << ", succes = " << count << ", len = " << lenbb << ", avg_iter = " << iteravg << endl;
		cout << ss.str();
	}
	void max_tabu(int ii) {
		//~ Benchmark::reset();
		qmax.push_back(vector<int>());
		int i, j, k, l, m, lbest;
		lbest = 0;
		int lenbest = 0;
		Titer = 0;
		int M_iter = 0;
		//~ starting_time = (double)clock();
		for( i = 0; i < len_time; i++ ) {
			l = tabu(len_improve);
			M_iter = M_iter + Iter; 
			if( l > lbest ) {
				lbest = l; 
				//~ finishing_time = (double)clock(); 
				Titer = M_iter; 
				len_W = len_best;       


				//~ dbgmsg("max.weight clique vertices = ");
				qmax.back().clear();
				//~ for (i = 0; i < len; i++) {
					//~ qmax.back().push_back(cruset[i]);
				//~ }
		        for( i = 0; i < Max_Vtx; i++ ) {
					if (vectex[i] > 0) {
						qmax.back().push_back(i);
						//~ dbgmsg("vectex[" << i << "] = " << vectex[i]);
					}
		            Tbest[ i ] = vectex[ i ];
		        }
		        verify();
		        /*for( i = 0; i < Max_Vtx; i++ )
		        {
		            Tbest[ i ] = vectex[ i ];
		        }*/
			}
			if( l == Waim ) goto save_lbest;
			//~ if( l == Waim ) return lbest;
		}
		save_lbest:
		W_used[ ii ] = lbest;
		len_used[ ii ] = len_W;
		Iteration[ ii ] = Titer;
		dbgmsg("i = " << ii << " l = " << lbest);
		//~ return lbest;
		//~ cout << "time to find one max.weight clique " << Benchmark::seconds_from_start() << " wallclock seconds" << endl;
	}
	void mnts(const int ii) {
		for(int i = 0; i < ii; i++ ) {
			max_tabu(i);
		}
	}
};


int main(int argc, char **argv) {
	try {
		if ( argc != 4 ) throw Error("die : usage ./progname filename(input graph file) Waim(objective weight value) len_improve(the search depth value)");
		MNTS mnts(argv[1], atoi(argv[2]), atoi(argv[3]));
		mnts.mnts(100);
		mnts.output();
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
}
//~ #endif
