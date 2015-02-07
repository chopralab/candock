#ifndef _MOTIF_H
#define _MOTIF_H

#include "const.h"


class Motif {
 private:
  struct cmp {
    bool operator()(string a, string b) {  // sortiramo stringe kot inte !
      return atoi(a.c_str()) < atoi(b.c_str());
    }
  };
  

  string delete_all_whitespaces(string);
  pair<size_t, size_t> ffo(string, set<string>, size_t=0);
  set<string> divide(string, string);
  set<ChResi> sele;
 public:
  static string debracket(char[]);
  Motif(string);
  bool is_in_motif(ChResi);
};

extern Motif *motif;

#endif // _MOTIF_H
