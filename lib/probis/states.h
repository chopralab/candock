/* MIT License
*
* Copyright (c) 2017 Janez Konc
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

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
