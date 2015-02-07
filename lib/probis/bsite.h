#ifndef _BSITE_H
#define _BSITE_H

#include "const.h"
#include "geo.h"
#include "ligand.h"

class Ligand;
class Coor;
class BSite;

bool by_bstype_bsscore(BSite*, BSite*);
bool by_bstype_bsscore_size(BSite*, BSite*);
bool by_bstype_size(BSite*, BSite*);


class BSite {  // opisuje en binding site
 public:
 BSite(Coor BSC, bstype BST) :  center(BSC), type(BST), score(0) {}

  Coor center;
  bstype type;
  float score;  // najboljsi cluster_score ?
  set<Ligand*, ligcomp_bsite> lig;
  set<ChResi> unija; // unija vseh aminokislin vseh ligandov za ta binding site
};




#endif // _BSITE_H
