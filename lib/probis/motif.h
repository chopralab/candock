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
