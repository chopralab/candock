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

#ifndef _RESIDUE_H
#define _RESIDUE_H

#include "const.h"

class ClusterData;
class Molecule;


class Residue {
 public:
// Residue(char CH, char RN, int RI, ClusterData *CD, bool FLX) : chain_id(CH), sig_id(0), resn(RN), resi(RI), cons(0), cd(CD), m(NULL), flx(FLX) {}
// Residue(char CH, char RN, int RI) : chain_id(CH), sig_id(0), resn(RN), resi(RI), cons(0), cd(NULL), m(NULL), flx(false) {}
// Residue(char CH, int RI) : chain_id(CH), sig_id(0), resn(' '), resi(RI), cons(0), cd(NULL), m(NULL), flx(false) {}
 Residue(char CH, char RN, int RI, ClusterData *CD, Molecule *M, bool FLX) : chain_id(CH), sig_id(0), resn(RN), resi(RI), cons(0), cd(CD), m(M), flx(FLX) {}
 Residue(char CH, char RN, int RI) : chain_id(CH), sig_id(0), resn(RN), resi(RI), cons(0), cd(NULL), m(NULL), flx(false) {}
 Residue(char CH, int RI) : chain_id(CH), sig_id(0), resn(' '), resi(RI), cons(0), cd(NULL), m(NULL), flx(false) {}
 Residue(char CH, int RI, Molecule *M) : chain_id(CH), sig_id(0), resn(' '), resi(RI), cons(0), cd(NULL), m(M), flx(false) {}
  
  char chain_id;
  int sig_id; // signature (fingerprint) residue number
  char resn;
  int resi;   // stevilka aminokisline v prileganem proteinu
  size_t rawcons, cons;
  ClusterData *cd;
  Molecule *m;
  bool flx;  // fleksibilna (a vseeno ohranjena v zaporedju) aminokislina
};


/* class operator za primerjavo Residue */
struct rescomp {
  bool operator() (const Residue*, const Residue*) const;
};



#endif // _RESIDUE_H
