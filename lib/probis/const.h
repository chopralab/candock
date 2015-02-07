#ifndef _CONST_H
#define _CONST_H

#include <stdio.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <mysql.h>

#include <thread>
#include <mutex>



#define XS      1e-60
#define SMALL   100
#define MEDIUM  1000
#define LARGE   10000
#define XLARGE  100000
#define XXL     100000000
#define XXXL    1e+20

typedef pair<char, int> ChResi;
typedef pair<ChResi,ChResi> TwoChResi;

enum state_type {TOTAL   = 0, 
//                 PPI     = 1, 
//                 BIND    = 2, 
                 WRITE   = 3, 
                 SURFDB  = 4, 
                 MARK    = 5, 
                 LIGAND  = 6, 
                 RESULTS = 7,
//                 SEQ     = 8,
                 ALIGN   = 9,
//                 BIOU    = 10,
                 HEAD    = 11,
                 ALLBIO  = 12
#ifdef CILE
                 ,PATCH  = 13
#endif
                 }; // BIND is not yet supported


enum {_firstModel = true, _allModels = false, _allChains = true, _selectedChains = false}; // argumenti za read_PDB
enum { _noHetero = 0, _hetero = 1, _nucleic = 2 };
enum bstype {_pp = 0, _nu = 1, _sl = 2, _io = 3}; // tip binding site-a
//enum { _noSaveMem = 0, _saveMemLig = 1, _saveMemBiou = 2};
enum { _noSaveMem = 0, _saveMemLig = 1};

/* za preveriti tip residueja */
const string amino = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL";
const string nucleic = "  A  T  G  U DA DT DG DC";
const string ions = " LI BE NA MG AL  K CA CR MN FE CO NI CU ZN PD PB AG CD SN SB CS BA PT AU HG UNK BR CL  F IOD";

/* napake */
class Err {
  const string name;
  const int code;
 public:
 Err(const string N, const int c) : name(N), code(c) {}
  string what() { return name; }
};

/* posebni znaki za json */
const string json_special = "\"\\";

/* PROGRAM OPTIONS & MODIFIERS */
extern state_type _state;
extern bool _noclus, _noprune, _srf, _verbose, _local, _nobb, _nomarkbb, _super, _caonly, _motif1, _motif2, _database, _pairwise, _nofp, _lig, _bsite1, _bsite2, _longnames;

/* DEFAULTS THAT CAN BE (SOME) OVERRIDEN IN PARAM FILE */
/* GENERAL VARIABLES */
extern int			NCPU				;
extern string      SRF_FILE               ;
extern string      PROTEIN1               ;
extern string      PROTEIN2               ;
extern string      CHAIN1                 ;   // chain1 and chain2 are both on protein1 for PPI, for other they are on prot1&2
extern string      CHAIN2                 ;
extern float       INTER_CHAIN_DIST       ;
//~ extern char        SURF_FILE[SMALL]       ;
extern string      SURF_FILE              ;
extern string      INDIR                  ;  // direktorij za input datoteke (razen tistih v zvezi z ligandi)
extern string      OUTDIR                 ;
/* SURF CONSTANTS */
#ifdef CILE
extern float       CILER                  ;
#endif
extern float       SPACING                ;
extern float       PROBE                  ;
extern float       MAXR                   ;
extern string      MOTIF1                 ;
extern string      MOTIF2                 ;
extern string      BSITE1                 ;
extern string      BSITE2                 ;
extern float       CATOM                  ;
extern float       NATOM                  ;
extern float       OATOM                  ;
extern float       SATOM                  ;
extern float       HATOM                  ;
extern float       EPS                    ;
extern float       SURF                   ;
extern float       DESC_PROBE_DIST        ;
extern float       BFACTOR                ;
extern float       POP_MEAN               ;   
extern float       POP_SD                 ;
extern float       NFP                    ;
//extern int         SIG_M                  ;   
//extern int         SIG_N                  ;
//extern float       SIG_CONS               ;
#define            SITE_RMSD        3.0
/* PRODUCT CONSTANTS */         
extern float       RESOLUTION             ;
extern float       OGA                    ;
extern float       OGB                    ;
extern float       OGC                    ;
extern float       OGV                    ;
extern float       THRMSD                 ;
/* DESC CONSTANTS */            
extern float       CUTOFF_FIRST           ;
extern float       CUTOFF_SECOND          ;
//extern float       PCUT                   ;
/* KABSCH CONSTANTS */          
extern float       NORM_EPS               ;
#define            MAX_VERTICES    20000
extern float       MATRIX_CUTOFF          ;
extern float       VECTOR_CUTOFF          ;
/* CLIQUE CONSTANTS */          
extern float       RMSD_INCR              ;
extern float       RMSD_INCR_ADD          ;
#define            WORD            32
extern float       CLIQUE_TLIMIT          ;          // time limit for clique calculation
/* CLUSTER CONSTANTS */                     
extern int         CLUS_SPEC              ;     
/* MYSQL CONSTANTS */           
//extern char        SERVER[15]             ;
//extern bool        PAIRWISE               ;             // not for mysql
//extern char        USER[15]               ;
//extern char        PASSWORD[15]           ;
//extern char        DATABASE[20]           ;
//extern int         WEIGHT                 ;     // OBSOLETE !!
extern int         CLUS_TO_OUTPUT         ;             // not for mysql
extern int         SCONS                  ;                // not for mysql
extern int         ALIGNMENT_NO                  ;                 // not for mysql
/* SCORE CONSTANTS */           
extern char        BLOSUM[15]             ;
extern float       SURF_VECTOR_ANGLE      ;         // pi/2
extern float       BLOSUM_SCORE           ;            // delete all with E-value more than this
extern float       CALPHA_RMSD            ;            // delete all with rmsd more than this
extern float       K                      ;
extern float       LAMBDA                 ;
extern float       Z_SCORE                ;
extern float       Z_SCORE_CONS           ;
extern float       Z_SCORE_FP             ;
/* GLOBAL ALIGNMENT CONSTANTS */
extern float       GLOBAL_DIST            ;
extern float       GLOBAL_SVAL            ;
extern float       GLOBAL_ANGLE           ;         // pi/4
/* LIGANDS CONSTANTS */
extern float       LIG_NEIGHB             ;
extern float       PP_DIST                ;
extern float       SL_DIST                ;
extern float       NU_DIST                ;
extern float       IO_DIST                ;
extern int         NCLTRIES               ;
extern int         NMOL                   ;
extern string      LIGDIR                 ; 
/* NOSQL CONSTANTS */
extern string      NOSQL_FILE             ; 
extern string      JSON_FILE             ; 
//extern string      NOSQL_DIR              ; 

///* RUSSELL TEST SET CONSTANTS */
//extern char        RUSSELL1[SMALL]        ;
//extern char        RUSSELL2[SMALL]        ;

extern string blosum80_3[25];
extern string blosum62[25];
extern string blosum80[25];

void close_resources();


#endif // _CONST_H
