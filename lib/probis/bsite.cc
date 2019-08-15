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
