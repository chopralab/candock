#ifndef _ARGS_H
#define _ARGS_H

#include "const.h"

class Args {
 public:
  Args();
  char *interface;      // not in use -- OCT/15/2008
  int t_cons_limit;     // not in use
  int dbsize;           // not in use 
  void fill_args(int, char*[]);
  void print_state();
  void print_modifiers();
  void print_constants();
  friend ostream& operator<< (ostream &out, Args &arg);
 private:
  map<string, string> opts;
  map<string, string> mods;
  void read_constants(string);
  void read_surf_inp_file(char*);
};

extern Args *args;



#endif // _ARGS_H
