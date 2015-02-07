#ifndef _BIT_H
#define _BIT_H

#include "const.h"


class Bit {
protected:
  bool get(int, int);
  void set_one(int, int);
  void set_zero(int, int);
  void set_zero_fast(int, int);
  unsigned int e[1+MAX_VERTICES/WORD][MAX_VERTICES];
};




#endif // _SCORE_H
