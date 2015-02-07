#ifndef _STATES_H
#define _STATES_H

#include "const.h"
#include "nosql.h"

class Args;


/* napovemo funkcije iz states.cc*/
void get_motif(Molecule*, string);
void get_bsite(Molecule*, string, string, string);
void state0(string, string);
void state0(NoSql::WriteStruct&, string, string);
//void state1(Args*);
//void state2(Args*);
void state3(Args*);
void state4(Args*);
void state5(Args*);
void state6(Args*);
void state7(Args*);
void state8(Args*);
void state9(Args*);
//void state10(Args*);
void state11(Args*);
void state12(Args*);
#ifdef CILE
void stateX(Args*);
#endif 

#endif // _STATES_H
