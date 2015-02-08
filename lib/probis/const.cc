#include "const.h"
#include "debug.hpp"
#include "args.h"
#include "geo.h"
#include "ligands.h"
#include "product.h"
#include "grid.h"
#include "probe.h"
#include "desc.h"
#include "clique.h"
#include "subgraph.h"
#include "cluster.h"
#include "motif.h"
#include "molecule.h"
#include "output.h"
#include "score.h"
#include "nosql.h"

Coor null_vect (0,0,0);

/* PROGRAM OPTIONS & MODIFIERS */
state_type _state;
bool _noclus=false, _noprune=false, _srf=false, _verbose=false, _local=false, _nobb=false, _nomarkbb=false, _super=false, _caonly=false, _motif1=false, _motif2=false, _database=false, _pairwise=false, _nofp=false, _lig=false, _bsite1=false, _bsite2=false, _longnames=false;

/* DEFAULTS THAT CAN BE (SOME) OVERRIDEN IN PARAM FILE */
/* GENERAL VARIABLES */
int			NCPU					= 1;
string      SRF_FILE               = "";
string      PROTEIN1               = "";
string      PROTEIN2               = "";
string      CHAIN1                 = "";   // chain1 and chain2 are both on protein1 for PPI, for other they are on prot1&2
string      CHAIN2                 = "";
float       INTER_CHAIN_DIST       = 3.0;
//~ char        SURF_FILE[SMALL]       = "";   // this is not written to mysql
string      SURF_FILE		       = "";   // this is not written to mysql
string      INDIR                  = "";
string      OUTDIR                 = "";
/* SURF CONSTANTS */
#ifdef CILE
float       CILER                  = 12.0;
#endif
float       SPACING                = 3.0;
float       PROBE                  = 1.4;
float       MAXR                   = 1.8;
string      MOTIF1                 = "";
string      MOTIF2                 = "";
string      BSITE1                 = "";
string      BSITE2                 = "";
//const int   NUM_CELLS              = 20;
//const int   NUM_CELLS_1            = 49; // NUM_CELLS - 1
float       CATOM                  = 1.5;
float       NATOM                  = 1.45;
float       OATOM                  = 1.35;
float       SATOM                  = 1.70;
float       HATOM                  = 1.0;
float       EPS                    = 0.001;
float       SURF                   = 1.1;
float       DESC_PROBE_DIST        = 3.0;
float       BFACTOR                = 0.8;
float       POP_MEAN               = 1.958;   
float       POP_SD                 = 2.21468;
float       NFP                    = 0.25;
//int         SIG_M                  = 1;   
//int         SIG_N                  = 5;
//float       SIG_CONS               = 7;
//map<int,int> SIG_RESI_CHECKED;
//int         SIG_N0                 = 5;
/* PRODUCT CONSTANTS */         
//float       MNSP_LOW               = 2.4;   // OBSOLETE
//float       MNSP_HIGH              = 3.3;   // OBSOLETE
float       RESOLUTION             = 2.0;
float       OGA                    = 1;
float       OGB                    = 1;
float       OGC                    = 1;
float       OGV                    = 2;
float       THRMSD                 = 1.2;
//float       RESOLUTION_2           = RESOLUTION*RESOLUTION;
//int         MAX_PRODUCT_GRAPH_SIZE = 1500000;
/* DESC CONSTANTS */            
//float   COULOMB                = 3.0;
float       CUTOFF_FIRST           = 9.0;
float       CUTOFF_SECOND          = 15.0;
//float       PCUT                   = 6.0;
/* KABSCH CONSTANTS */          
float       NORM_EPS               = 0.00000001;
//const int   MAX_VERTICES           = 20000;
float       MATRIX_CUTOFF          = 0.2;
float       VECTOR_CUTOFF          = 2.0;
/* CLIQUE CONSTANTS */          
float       RMSD_INCR              = 0.50;
float       RMSD_INCR_ADD          = 0.60;
//const int   WORD                   = 32;
float       CLIQUE_TLIMIT          = 1.0;          // time limit for clique calculation
//int         MIN_CLIQUE_SIZE        = 3;
//  /* PATCH CONSTANTS */           
//  float       MAX_INTER_DIST         = 8.0;   // OBSOLETE !!
//  int         MIN_DENSITY            = 2;     // OBSOLETE !!
//  const int   MAX_INTERFACE          = 1000;  // OBSOLETE !!
//  const int   MAX_NEIGHB             = 50;    // OBSOLETE !!
/* CLUSTER CONSTANTS */                     
int         CLUS_SPEC              = 5;     
/* MYSQL CONSTANTS */           
//char        SERVER[15]             = "localhost";
//bool        PAIRWISE               = false;             // not for mysql
//char        USER[15]               = "root";
//char        PASSWORD[15]           = "";
//char        DATABASE[20]           = "results_db";
//int         WEIGHT                 = 5;     // OBSOLETE !!
int         CLUS_TO_OUTPUT         = 5;             // not for mysql
int         SCONS                  = 4;                // not for mysql
int         ALIGNMENT_NO           = 0;                 // not for mysql
/* SCORE CONSTANTS */           
char        BLOSUM[15]             = "blosum80";
float       SURF_VECTOR_ANGLE      = 1.5708;         // pi/2
float       BLOSUM_SCORE           = 1.0e+4;            // delete all with E-value more than this
float       CALPHA_RMSD            = 3.0;            // delete all with rmsd more than this
float       K                      = 0.17;
float       LAMBDA                 = 0.343;
float       Z_SCORE                = 1.0;
float       Z_SCORE_CONS           = 2.0;
float       Z_SCORE_FP             = 3.0;
/* GLOBAL ALIGNMENT CONSTANTS */
float       GLOBAL_DIST            = 10.0;
float       GLOBAL_SVAL            = 23;
float       GLOBAL_ANGLE           = 0.7854;         // pi/4
/* LIGANDS CONSTANTS */
float       LIG_NEIGHB             = 2.00;
float       PP_DIST                = 10.00;
float       SL_DIST                = 8.00;
float       NU_DIST                = 10.00;
float       IO_DIST                = 5.00;
//float       PP_DIST                = 7.00;
//float       SL_DIST                = 5.00;
//float       NU_DIST                = 6.00;
//float       IO_DIST                = 5.00;
int         NCLTRIES               = 100; 
int         NMOL                   = 0;     // iz koliko prvih prileganih proteinov bo izracunalo ligande
string      LIGDIR                 = "";   // direktorij s prileganimi ligandi (.lig)
/* NOSQL CONSTANTS */
//string      NOSQL_DIR               = ".";   // direktorij z nosql bazo
string      NOSQL_FILE       = "";
string      JSON_FILE        = ""; 




string blosum80_3[25] = { 
  "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *", 
  "A  7 -3 -3 -3 -1 -2 -2  0 -3 -3 -3 -1 -2 -4 -1  2  0 -5 -4 -1 -3 -2 -1 -8",  
  "R -3  9 -1 -3 -6  1 -1 -4  0 -5 -4  3 -3 -5 -3 -2 -2 -5 -4 -4 -2  0 -2 -8",  
  "N -3 -1  9  2 -5  0 -1 -1  1 -6 -6  0 -4 -6 -4  1  0 -7 -4 -5  5 -1 -2 -8",  
  "D -3 -3  2 10 -7 -1  2 -3 -2 -7 -7 -2 -6 -6 -3 -1 -2 -8 -6 -6  6  1 -3 -8",  
  "C -1 -6 -5 -7 13 -5 -7 -6 -7 -2 -3 -6 -3 -4 -6 -2 -2 -5 -5 -2 -6 -7 -4 -8",  
  "Q -2  1  0 -1 -5  9  3 -4  1 -5 -4  2 -1 -5 -3 -1 -1 -4 -3 -4 -1  5 -2 -8",  
  "E -2 -1 -1  2 -7  3  8 -4  0 -6 -6  1 -4 -6 -2 -1 -2 -6 -5 -4  1  6 -2 -8",  
  "G  0 -4 -1 -3 -6 -4 -4  9 -4 -7 -7 -3 -5 -6 -5 -1 -3 -6 -6 -6 -2 -4 -3 -8",  
  "H -3  0  1 -2 -7  1  0 -4 12 -6 -5 -1 -4 -2 -4 -2 -3 -4  3 -5 -1  0 -2 -8",  
  "I -3 -5 -6 -7 -2 -5 -6 -7 -6  7  2 -5  2 -1 -5 -4 -2 -5 -3  4 -6 -6 -2 -8",  
  "L -3 -4 -6 -7 -3 -4 -6 -7 -5  2  6 -4  3  0 -5 -4 -3 -4 -2  1 -7 -5 -2 -8",  
  "K -1  3  0 -2 -6  2  1 -3 -1 -5 -4  8 -3 -5 -2 -1 -1 -6 -4 -4 -1  1 -2 -8",  
  "M -2 -3 -4 -6 -3 -1 -4 -5 -4  2  3 -3  9  0 -4 -3 -1 -3 -3  1 -5 -3 -2 -8",  
  "F -4 -5 -6 -6 -4 -5 -6 -6 -2 -1  0 -5  0 10 -6 -4 -4  0  4 -2 -6 -6 -3 -8",  
  "P -1 -3 -4 -3 -6 -3 -2 -5 -4 -5 -5 -2 -4 -6 12 -2 -3 -7 -6 -4 -4 -2 -3 -8",  
  "S  2 -2  1 -1 -2 -1 -1 -1 -2 -4 -4 -1 -3 -4 -2  7  2 -6 -3 -3  0 -1 -1 -8",  
  "T  0 -2  0 -2 -2 -1 -2 -3 -3 -2 -3 -1 -1 -4 -3  2  8 -5 -3  0 -1 -2 -1 -8",  
  "W -5 -5 -7 -8 -5 -4 -6 -6 -4 -5 -4 -6 -3  0 -7 -6 -5 16  3 -5 -8 -5 -5 -8",  
  "Y -4 -4 -4 -6 -5 -3 -5 -6  3 -3 -2 -4 -3  4 -6 -3 -3  3 11 -3 -5 -4 -3 -8",  
  "V -1 -4 -5 -6 -2 -4 -4 -6 -5  4  1 -4  1 -2 -4 -3  0 -5 -3  7 -6 -4 -2 -8",  
  "B -3 -2  5  6 -6 -1  1 -2 -1 -6 -7 -1 -5 -6 -4  0 -1 -8 -5 -6  6  0 -3 -8",  
  "Z -2  0 -1  1 -7  5  6 -4  0 -6 -5  1 -3 -6 -2 -1 -2 -5 -4 -4  0  6 -1 -8",  
  "X -1 -2 -2 -3 -4 -2 -2 -3 -2 -2 -2 -2 -2 -3 -3 -1 -1 -5 -3 -2 -3 -1 -2 -8",  
  "* -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1"  
};


string blosum62[25] = { 
  "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *", 
  "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4",  
  "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4",  
  "N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4",  
  "D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4",  
  "C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4",  
  "Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4",  
  "E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4",  
  "G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4",  
  "H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4",  
  "I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4",  
  "L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4",  
  "K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4",  
  "M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4",  
  "F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4",  
  "P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4",  
  "S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4",  
  "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4",  
  "W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4",  
  "Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4",  
  "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4",  
  "B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4",  
  "Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4",  
  "X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4",  
  "* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1"
}; 


string blosum80[25] = { 
  "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *",   
  "A  5 -2 -2 -2 -1 -1 -1  0 -2 -2 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -6",  
  "R -2  6 -1 -2 -4  1 -1 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -4 -3 -3 -2  0 -1 -6",  
  "N -2 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4  4  0 -1 -6",  
  "D -2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -6 -4 -4  4  1 -1 -6",  
  "C -1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1 -4 -4 -1 -6",  
  "Q -1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3  0  3 -1 -6",  
  "E -1 -1 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3  1  4 -1 -6",  
  "G  0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4 -1 -3 -1 -6",  
  "H -2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  2 -4 -1  0 -1 -6",  
  "I -2 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3 -4 -4 -1 -6",  
  "L -2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -3 -3 -2 -2 -2  1 -4 -3 -1 -6",  
  "K -1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3 -1  1 -1 -6",  
  "M -1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1 -3 -2 -1 -6",  
  "F -3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1 -4 -4 -1 -6",  
  "P -1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -3 -1 -3 -4  8 -1 -2 -5 -4 -3 -2 -2 -1 -6",  
  "S  1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2  0  0 -1 -6",  
  "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0 -1 -1 -1 -6",  
  "W -3 -4 -4 -6 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -5 -4 -4 11  2 -3 -5 -4 -1 -6",  
  "Y -2 -3 -3 -4 -3 -2 -3 -4  2 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2 -3 -3 -1 -6",  
  "V  0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  4 -4 -3 -1 -6",  
  "B -2 -2  4  4 -4  0  1 -1 -1 -4 -4 -1 -3 -4 -2  0 -1 -5 -3 -4  4  0 -1 -6",  
  "Z -1  0  0  1 -4  3  4 -3  0 -4 -3  1 -2 -4 -2  0 -1 -4 -3 -3  0  4 -1 -6",  
  "X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -6",  
  "* -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1"
};

