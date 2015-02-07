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
