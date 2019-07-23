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

#ifndef _GEO_H
#define _GEO_H

#include "const.h"

class Atom;

template <class T>

inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}



class Coor {
 public:
  Coor() {}
  Coor (double X, double Y, double Z) : x(X), y(Y), z(Z) {}
  double x, y, z;
};

float z_score(float);
vector<string> reg_split(string, string, string="", string="");
string json_escape_special_characters(string);
string add_end_slash(string);
string remove_whitespaces(string);
float str2float(const string&);
double ncomb(int, int);
double significance(int, int, int, int);
void to_upper (char*);
//void extract_pdb_code(char*);
string extract_pdb_code(const string);
char one_letter_code(char*);
inline double dist(Coor v1, Coor v2) {   return sqrt(pow(v1.x-v2.x,2) + pow(v1.y-v2.y,2) + pow(v1.z-v2.z,2)); }
inline double dist_fast(Coor v1, Coor v2) {  return pow(v1.x-v2.x,2) + pow(v1.y-v2.y,2) + pow(v1.z-v2.z,2); }
Coor operator% (Coor, Coor);
Coor operator* (Coor, double);
Coor operator- (Coor, double);
Coor operator+ (Coor, double);
Coor operator- (Coor, Coor);
Coor operator+ (Coor, Coor);
double operator* (Coor, Coor);
bool operator== (Coor, Coor);
Coor project(Coor, Coor);
void set_zero(Coor&);
Coor norm(Coor);
Coor intersection(Coor, Coor, Coor, Coor);
double refpoint(double, double, double, double);   // ta funkcija je definirana kot 'inline'
bool probe_center(Atom*, Atom*, Atom*, Coor&, Coor&);
double angle(Coor, Coor, Coor);
Coor rotate_vector(Coor, Coor, Coor, double, double);
Coor calculate_f2(Atom, Coor, double);
double line_dist(Coor, Coor, Coor, Coor);
double line_point_dist(Coor, Coor, Coor);
bool check_direction(Atom, Atom, Atom, Coor);


extern Coor null_vect;

#endif  // _GEO_H
