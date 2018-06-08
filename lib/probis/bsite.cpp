#include "bsite.h"

bool by_bstype_bsscore(BSite *i, BSite *j) {
  //  type  bs-score
  //  ------------------------------------
  //  _pp     7.1              _pp     7.1   --> best of type _pp
  //  _ion    7.5     --->     _pp     6.0
  //  _pp     6.0              _ion    7.5   --> best of type _ion
  //  _ion    3.5              _ion    3.5
  //
  if (i->type < j->type) 
    return true;
  else if ( i->type == j->type) 
    if (i->score > j->score) 
      return true;
  return false;

}

bool by_bstype_bsscore_size(BSite *i, BSite *j) {
  /* 
     najprej posortira skupaj enak tipe binding site-ov, nato sortira po
     bs-score in znotraj enakih bs-score po stevilu vseh aminokislin v tem
     binding site-u
   */
  if (i->type < j->type) {
    return true;
  }
  else if ( i->type == j->type) {
    if (i->score > j->score) {
      return true;
    }
    else if ( i->score == j->score) {
      if (i->unija.size() > j->unija.size())
        return true;
    }
  }
  return false;

}

bool by_bstype_size(BSite *i, BSite *j) {
  /* 
     Najprej posortira skupaj enak tipe binding site-ov, nato sortira po
     po stevilu vseh aminokislin v tem binding site-u.
  */
  if (i->type < j->type) {
    return true;
  }
  else if ( i->type == j->type) {
    if (i->unija.size() > j->unija.size())
      return true;
  }
  return false;
}
