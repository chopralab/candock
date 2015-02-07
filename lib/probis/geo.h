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
