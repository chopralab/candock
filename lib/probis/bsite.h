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
