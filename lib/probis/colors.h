#ifndef _COLORS_H
#define _COLORS_H

class Colors {
 public:
 Colors() : size(0) {}
  int color[10];
  int size;
  bool is_visited(int);
  void copy(Colors);
};

#endif // _COLORS_H
