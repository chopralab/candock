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

#ifndef _BIOUNIT_H
#define _BIOUNIT_H


class BioUnit {
 public:
  ~BioUnit() {
    /* sprostimo pomnilnik za rotacijske matrike in translacijske vektorje biounitov */
    for (map<int, gsl_matrix*>::iterator it = U.begin(); it != U.end(); it++) {
      gsl_matrix_free(it->second);
    }
    for (map<int, gsl_vector*>::iterator it = t.begin(); it != t.end(); it++)
      gsl_vector_free(it->second);
  }

  map<int, gsl_matrix*> U;
  map<int, gsl_vector*> t;
//  set<char> chain_id;
  multimap<int, char> chain_id;
  vector<string> rota_pdb; // pdb file v tekstovni obliki
  string filename;

  bool najdi(string c_id) {
//    for (set<char>::iterator it = chain_id.begin(); it != chain_id.end(); it++) 
    for (multimap<int,char>::iterator it = chain_id.begin(); it != chain_id.end(); it++) 
      if (c_id.find(it->second) != string::npos) 
        return true; 
    return false; 
  }

  bool is_in_rota_num(int rota_num, char c_id) {
    pair<multimap<int,char>::iterator,multimap<int,char>::iterator> ret = chain_id.equal_range(rota_num);
    for (multimap<int,char>::iterator it = ret.first; it != ret.second; it++) 
      if (c_id == it->second) 
        return true; 
    return false; 
  }

};




#endif // _BIOUNIT_H
