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

#include "item.h"
#include "product.h"


void Item::copy_unique(Item *item) {

  for (unsigned int i = 0; i < item->vert.size(); i++) {
    unsigned int j;
    for (j = 0; j < vert.size(); j++)
//      if (vert[ j ]->row == item->vert[ i ]->row)
      if ((vert[ j ]->row->desc1 == item->vert[ i ]->row->desc1) && (vert[ j ]->row->desc2 == item->vert[ i ]->row->desc2))
        break;
    if (j == vert.size() ) 
      vert.push_back(new Vert(0, 0, 0, item->vert[ i ]->row ) );
//      vert[size++].row = item->vert[i].row;
  }
}

