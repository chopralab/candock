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

