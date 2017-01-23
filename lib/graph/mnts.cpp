#include "mnts.hpp"
void MNTS::initialize() {
	srand( (unsigned) time( NULL ) );
	Iteration = unique_ptr<int[]>(new int [iter]);
	len_used = unique_ptr<int[]>(new int [iter]);
	W_used = unique_ptr<int[]>(new int [iter]);
	vectex = unique_ptr<int[]>(new int [Max_Vtx]);
	funch  = unique_ptr<int[]>(new int [Max_Vtx]);
	address= unique_ptr<int[]>(new int [Max_Vtx]);
	tabuin = unique_ptr<int[]>(new int [Max_Vtx]);
	C0 = unique_ptr<int[]>(new int[Max_Vtx]);
	C1 = unique_ptr<int[]>(new int[Max_Vtx]);
	TC1= unique_ptr<int[]>(new int[Max_Vtx]);
	BC = unique_ptr<int[]>(new int[Max_Vtx]);
	FC1= unique_ptr<int[]>(new int[Max_Vtx]);
	cruset = unique_ptr<int[]>(new int [2000]);
	//~ Tbest = unique_ptr<int[]>(new int[Max_Vtx];
	//~ TTbest = unique_ptr<int[]>(new int[Max_Vtx]);
	for( int x = 0 ; x < Max_Vtx ; x++ ) {
		vectex[x]  = x;
		address[x] = x;
		BC[ x ] = 0;
		//~ Edge[x][x] = false;
		Edge[x][x] = true;
	}
	
	adjaclen = unique_ptr<int[]>(new int[Max_Vtx]);
	adjacMatrix = unique_ptr<unique_ptr<int[]>[]>(new unique_ptr<int[]>[Max_Vtx]);
	dbgmsg("max_weight_clique --- start filling : Size of graph = " << Max_Vtx);
	for( int x = 0 ; x < Max_Vtx ; x++ ) {
		adjaclen[x] = 0;
		adjacMatrix[x] = unique_ptr<int[]>(new int[Max_Vtx]);
		for( int y = 0 ; y < Max_Vtx ; y++ ) {
			if (!Edge[x][y]) {
				adjacMatrix[x][adjaclen[x]] = y;
				//~ dbgmsg("max_weight_clique : x = " << x << " y = " << y
					//~ << " adjacMatrix = " << adjacMatrix[x][adjaclen[x]]);

				adjaclen[x]++;
			}
		}
	}
	dbgmsg("max_weight_clique --- end filling");
	

}
void MNTS::clearGamma() {
	int tm1 = Max_Vtx*sizeof( int );
	memset( vectex.get(), 0, tm1 );
	memset(  funch.get(), 0, tm1 );
	memset(address.get(), 0, tm1 );
	memset( tabuin.get(), 0, tm1 );
	for( int i = 0; i < Max_Vtx; i++ ) {
		C0[ i ] = i;
		address[ i ] = i;
	}
	len0 = Max_Vtx;
	len1 = 0;
	len = 0;
	Wf = 0;
	Wbest = 0;
}
int MNTS::selectC0() {
	int i, k, l;
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
int MNTS::WselectC0( ) {
	int i, k, l1, l2, w1, w2;
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
int MNTS::expand(int SelN)
{
    int i, k1, m, n, n1;
    
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
        /*for( i = 0; i < Max_Vtx; i++ )
        {
            Tbest[ i ] = vectex[ i ];
        }*/
     }
    
    return 1;   
}

int MNTS::selectC1( )
{
    int i, k, l;
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

int MNTS::WselectC1( )
{
     int i, j, k=0, l, l1, l2, wmn, w1, w2, m, n;
     l1 = 0;
     l2 = 0;
     w1 = -1000000;
     w2 = -1000000;
     l = 0;
     for( i = 0; i < len1; i++ )
     {
         m = C1[ i ];
         n = BC[ m ];
         //~ if( (vectex[ n ] == 1) && (Edge[ m ][ n ] == 1) )
         //~ if( (vectex[ n ] == 1) && Edge[ m ][ n ] )
         if( (vectex[ n ] == 1) && !Edge[ m ][ n ] )
           l++;
         else
         {
             for( j = 0; j < len; j++ )
             {
                k = cruset[ j ];
                if( !Edge[ m ][ k ] )
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

int MNTS::plateau( int SelN )
{
     int i, k1, m, m1=0, n, n1, ti;
     
     m = C1[ SelN  ];
     for(ti = 0; ti < len; ti++)
     {
         m1 = cruset[ ti ];
         //~ if( Edge[ m1 ][ m ] == 1 )
         //~ if( Edge[ m1 ][ m ] )
         if( !Edge[ m1 ][ m ] )
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
        /*for( i = 0; i < Max_Vtx; i++ )
        {
            Tbest[ i ] = vectex[ i ];
        }*/
     }
     return 1;   
}


int MNTS::Mumi_Weigt()
{
    int i, k, l1;
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

int MNTS::backtract()
{
     int i, m1, n, ti, k1, n1;
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
     
     return 0;
}

int MNTS::tabu( int Max_Iter )
{
     int k, am, am1, ww, ww1, ww2, ti, m1;
     Iter = 0;
     clearGamma(); 
     while( 1 )
     {
        am = selectC0();
        if( am != -1 )
        {
            expand( am );
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
                expand( am );
                
                Iter++;
                if( Wbest == Waim )
                   return Wbest;
            }
            else
            {
                plateau( am1 );
                if( Wbest == Waim )
                    return Wbest; 
                Iter++;
            }
        }
        else if( (am != -1) && (am1 == -1) )
        {
             expand( am );
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
                plateau( am1 );
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
//~ void verify()
//~ {
     //~ int i, j, k1, k2, l, m;
     //~ for( i = 0; i < Max_Vtx; i++ )
     //~ {
          //~ if( Tbest[ i ] == 1 )
          //~ {
              //~ for( j = i+1; j < Max_Vtx; j++ )
              //~ if( (Tbest[ j ] == 1) && !Edge[ i ][ j ] )
                  //~ throw Error("hello there is something wrong");
          //~ }
     //~ }
//~ }
void MNTS::output() {
	int i; 
	stringstream ss;
	//~ for( i = 0; i < 100; i++ ) {
	for( i = 0; i < iter; i++ ) {
		ss << "sum = " << W_used[ i ] << ", iter = " << Iteration[ i ] << ", len = " << len_used[ i ] << ", qmax = "; 
		for (auto &v : qmax[i]) ss << v << " ";
		ss << endl;
	}
	ss << "\n\n the total information:" << endl;
	int wavg, iteravg, lenbb, success;
	wavg = iteravg = lenbb = success = 0;
	int best_v = 0;
	//~ for( i = 0; i < 100; i++ ) {
	for( i = 0; i < iter; i++ ) {
		if( W_used[ i ] > best_v ) {
			best_v = W_used[ i ];  
			lenbb = len_used[ i ];
		}
    }
	int count = 0;
	ss << "\n The best weight value for the maximum weighted problem is " << best_v << endl;
	//~ for( i = 0; i < 100; i++ ) {
	for( i = 0; i < iter; i++ ) {
		wavg = wavg + W_used[ i ];
	}  
	//~ double twavg = (double (wavg)) / 100 ; 
	double twavg = (double (wavg)) / iter ; 
	//~ for( i = 0; i < 100; i++ ) {
	for( i = 0; i < iter; i++ ) {
		if( W_used[ i ] == best_v ) {
			count++;
			iteravg = iteravg + Iteration[ i ];
		}
	}
	iteravg =  int ( (double (iteravg)) / count );
	ss << "avg_sum = " << twavg << ", succes = " << count << ", len = " << lenbb << ", avg_iter = " << iteravg << endl;
	dbgmsg(ss.str());
}
void MNTS::max_tabu(int ii) {
	//~ Benchmark::reset();
	qmax.push_back(vector<int>());
	int i, l, lbest;
	lbest = 0;
	Titer = 0;
	int M_iter = 0;
	//~ starting_time = (double)clock();
	for( i = 0; i < len_time; i++ ) {
		l = tabu(len_improve);
		//~ dbgmsg("i = " << i << " l = " << l << " lbest = " << lbest << " Waim = " << Waim << " len_improve = " << len_improve << " len_time = " << len_time);
		M_iter = M_iter + Iter; 
		if( l > lbest ) {
			lbest = l; 
			//~ finishing_time = (double)clock(); 
			Titer = M_iter; 
			len_W = len_best;       
			qmax.back().clear();
	        for( i = 0; i < Max_Vtx; i++ ) {
				if (vectex[i] > 0) {
					qmax.back().push_back(i);
					//~ dbgmsg("vectex[" << i << "] = " << vectex[i]);
				}
	            //~ Tbest[ i ] = vectex[ i ];
	        }
	        //~ verify();
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
