#include "geo.h"
#include "atom.h"
#include <cmath>

float z_score(float raw_score) {
  return (raw_score - POP_MEAN) / POP_SD;
}

/*
 * Helper function for splitting string to string chunks. Delimiter je lahko string, npr "),(",
 * dolocena morata biti prvi in zadnji znak v stringu s (i.e., f in l):
 * primer: (1,2,3),(2,4,6),(8,10,12) 
 * delim="),("
 * f="("
 * l=")"
 */
vector<string> reg_split(string s, string delim, string f, string l) {
  vector<string> ret_vector;
  
  if (f.length() < s.length()) {
    size_t last_i = f.length();
    size_t i = 0;
    while (i != string::npos) {
      i = s.find(delim, last_i);
      string part = (i==string::npos ? s.substr(last_i, s.length() - l.length() - last_i) : s.substr(last_i, i - last_i));
      ret_vector.push_back(part);
      last_i = i + delim.length();
    }
  }
  return ret_vector;
}

string json_escape_special_characters(string str) {
  /*
    V str-ju eskejpamo posebne znake, ki jih json parser ne prepozna
  */
  size_t found=str.find_first_of(json_special);
  while (found!=string::npos) {
        str.insert(found, "\\");
        //        cout << "Warning (Ligands) : Found special character [" << str.at(found) << "] in json string : " << str << endl; 
        //        str.erase(found, 1);
//        cout << "Warning (Ligands) : Found special character [" << str.at(found+1) << "] in json string : " << str << endl; 
        found=str.find_first_of(json_special, found+2);
      }
  return str;
}

string add_end_slash(string str) {
  /*
    Dodamo (ce ga se ni) zakljucni 'slash' simbol v ime direktorija.
  */
  if (str.size() > 0) {
    size_t found=str.find_last_of("/");
    if (found != str.size() - 1) str += "/";
  }
  return str;
}

string remove_whitespaces(string str) {
  /* 
     Zbrisemo blank znake na koncu in na zacetku stringa.
  */
  string whitespaces (" \t\f\v\n\r");
  size_t found;
  // od zadaj
  found=str.find_last_not_of(whitespaces);
  if (found!=string::npos)
    str.erase(found+1);
  else
    str.clear();            // str je cel whitespace
  // od spredaj
  found=str.find_first_not_of(whitespaces);
  if (found!=string::npos)
    str.erase(0, found);
  else
    str.clear();            // str je cel whitespace

  return str;

}


float str2float(const string& str) {
  stringstream ss;
  ss << str;
  float f;
  ss >> f;
  return f;
}

double ncomb(int n, int k) {
  double nk = 1;
  for (int i = 0; i < k; i++) {
    nk = nk * (n - i) / (k - i);
  }
  return nk;
}

double significance(int Ss, int Ps, int Is, int Os) {
  double sig = 0; 
  int min;
  if (Ps < Is) min = Ps; else min = Is;

  for (int i = Os; i < min; i++) {
    sig = sig + ncomb(Is, i) * ncomb(Ss - Is, Ps - i);
  }

  sig = sig / ncomb(Ss, Ps);

  return sig;
}

void to_upper (char *str) {
  int i=0;
  char c;
  while (str[i]) {
    c=str[i];
    str[i] = toupper(c);
    i++;
  }
}

string extract_pdb_code(const string fn) {
  /*
    Extract pdb code from a filename, e.g., /home/konc/1aze.pdb will output 1aze.
    Note, that previous version also converted pdb code to uppercase. Now this 
    does not happen and the pdb code is case sensitive!
  */

  size_t found1, found2;
  found1 = fn.find_last_of("/\\");  // dela tudi pod windowsi :)
  found2 = fn.find_last_of(".");  // najdemo suffix
  return fn.substr(found1 + 1, found2 - found1 - 1);
}

char one_letter_code(char *code) {
  to_upper(code);
  if (!strcmp(code, "ALA")) return 'A';
  if (!strcmp(code, "CYS")) return 'C';
  if (!strcmp(code, "ASP")) return 'D';
  if (!strcmp(code, "GLU")) return 'E';
  if (!strcmp(code, "PHE")) return 'F';
  if (!strcmp(code, "GLY")) return 'G';
  if (!strcmp(code, "HIS")) return 'H';
  if (!strcmp(code, "ILE")) return 'I';
  if (!strcmp(code, "LYS")) return 'K';
  if (!strcmp(code, "LEU")) return 'L';
  if (!strcmp(code, "MET")) return 'M';
  if (!strcmp(code, "ASN")) return 'N';
  if (!strcmp(code, "PRO")) return 'P';
  if (!strcmp(code, "GLN")) return 'Q';
  if (!strcmp(code, "ARG")) return 'R';
  if (!strcmp(code, "SER")) return 'S';
  if (!strcmp(code, "THR")) return 'T';
  if (!strcmp(code, "VAL")) return 'V';
  if (!strcmp(code, "TRP")) return 'W';
  if (!strcmp(code, "TYR")) return 'Y';
  return 'X';
}

Coor operator% (Coor v1, Coor v2) { // vektorski produkt    
  Coor res;
  res.x = v1.y*v2.z - v1.z*v2.y;
  res.y = v1.z*v2.x - v1.x*v2.z;
  res.z = v1.x*v2.y - v1.y*v2.x;
  return res;
}

Coor operator* (Coor v, double r) {   // mnozenje vektorja s konstanto
  Coor res;
  res.x = v.x*r;
  res.y = v.y*r;
  res.z = v.z*r;
  return res;
}

Coor operator- (Coor v1, double r) { // odstej konstano od vektorja
  Coor res;
  res.x = v1.x - r;
  res.y = v1.y - r;
  res.z = v1.z - r;
  return res;
}

Coor operator+ (Coor v1, double r) { // pristej konstano k vektorju
  Coor res;
  res.x = v1.x + r;
  res.y = v1.y + r;
  res.z = v1.z + r;
  return res;
}


Coor operator- (Coor v1, Coor v2) { // odstevanje vektorjev
  Coor res;
  res.x = v1.x - v2.x;
  res.y = v1.y - v2.y;
  res.z = v1.z - v2.z;
  return res;
}

Coor operator+ (Coor v1, Coor v2) { // sestevanje vektorjev
  Coor res;
  res.x = v1.x + v2.x;
  res.y = v1.y + v2.y;
  res.z = v1.z + v2.z;
  return res;
}

double operator* (Coor v1, Coor v2) {  // skalarni produkt
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

bool operator== (Coor v1, Coor v2) {
  return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

Coor project(Coor v1, Coor v2) { // project v1 on v2
  double dv2 = v2 * v2;
  double dv1v2 = v1 * v2;
  Coor v3;
  v3.x = v2.x * dv1v2 / dv2;
  v3.y = v2.y * dv1v2 / dv2;
  v3.z = v2.z * dv1v2 / dv2;
  return v3;
}

void set_zero(Coor &v) {
  v.x = 0.0;
  v.y = 0.0;
  v.z = 0.0;
}


Coor norm(Coor v) {
  Coor res;
  double d = sqrt(pow(v.x,2) + pow(v.y,2) + pow(v.z,2));
  if (d == 0) { 
    printf("zero length vector!\n"); 
    return null_vect; 
  }
  res.x = v.x/d;
  res.y = v.y/d;
  res.z = v.z/d;
  return res;
}

Coor intersection(Coor t1, Coor e, Coor t2, Coor f) { // presecisce dveh premic
  /*
    We calculate the intersection of two lines defined by vector point t1 and directional vector e
    and vector point t2 and directional vector f. We have to divide with the greatest tmp{1,2,3}.
  */
  double lambda1;
  double tmp1, tmp2, tmp3;
  tmp1 = e.z * f.x - e.x * f.z;
//  tmp1 = e.x * f.y - e.y * f.x;
  tmp2 = e.z * f.y - e.y * f.z;
  tmp3 = e.y * f.x - e.x * f.y;
  if (fabs(tmp1) > fabs(tmp2))
    if (fabs(tmp1) > fabs(tmp3))
        lambda1 = (f.z * (t1.x - t2.x) + f.x * (t2.z - t1.z)) / tmp1;
  //        lambda1 = (f.x * (t1.y - t2.y) + f.y * (t2.x - t1.x)) / tmp1;
    else
      lambda1 = (f.y * (t1.x - t2.x) + f.x * (t2.y - t1.y)) / tmp3;
  else
    if (fabs(tmp2) > fabs(tmp3))
      lambda1 = (f.z * (t1.y - t2.y) + f.y * (t2.z - t1.z)) / tmp2;
    else
      lambda1 = (f.y * (t1.x - t2.x) + f.x * (t2.y - t1.y)) / tmp3;

  return t1 + e * lambda1;
}

inline double refpoint(double r1, double r2, double r3, double dist12) {
  /*
    On the line between the centers of atoms at1 and at2 is point X which is directly below 
    the center of the probe sphere as it touches the two atoms. This function returns the 
    distance between point X and atom center at1.
  */
  
  return (pow(r1+r3,2) - pow(r2+r3,2) + pow(dist12,2)) / (2*dist12); 
}

bool probe_center(Atom *at1, Atom *at2, Atom *at3, Coor &center_1, Coor &center_2) {
/*
  In this function we calculate the spherical probe center as it is in between the 
  three atoms at1, at2, at3;
  Input are the three atoms.
  Output is the two center vectors (up and below atoms at1, at2, at3); the function
  returns false if the center of the probe atom cannot be calculated.
*/
  Coor n, X1, X2, c1_c2, c1_c3, p1, p2, tmp;
  Coor Y;
  double dT4Y;
  c1_c2 = at2->crd - at1->crd;
  c1_c3 = at3->crd - at1->crd;
  n = c1_c2 % c1_c3;
//tukaj ostal
//  at1->print_atom();
//  at2->print_atom();
//  at3->print_atom();
//  cout << "pred norm" << endl;
//  cout << "n = " << n.x << " " << n.y << " " << n.z << endl;
  n = norm(n);
//  cout << "n = " << n.x << " " << n.y << " " << n.z << endl;
//  cout << "c1_c2 = " << c1_c2.x << " " << c1_c2.y << " " << c1_c2.z << endl;
  c1_c2 = norm(c1_c2);
//  cout << "c1_c2 = " << c1_c2.x << " " << c1_c2.y << " " << c1_c2.z << endl;
//  cout << "c1_c3 = " << c1_c3.x << " " << c1_c3.y << " " << c1_c3.z << endl;
  c1_c3 = norm(c1_c3);  
//  cout << "c1_c3 = " << c1_c3.x << " " << c1_c3.y << " " << c1_c3.z << endl;
//  cout << "za norm" << endl;
  X1 = at1->crd + c1_c2 * refpoint(at1->r, at2->r, PROBE, dist(at1->crd, at2->crd));
  X2 = at1->crd + c1_c3 * refpoint(at1->r, at3->r, PROBE, dist(at1->crd, at3->crd));
  p1 = c1_c2 % n;
  p2 = c1_c3 % n;
  Y = intersection(X1 , p1, X2, p2);
  tmp = Y - at1->crd;
  dT4Y = pow(at1->r + PROBE, 2) - tmp * tmp;
  if (dT4Y < 0) return false;
  dT4Y = sqrt(dT4Y);
  center_1 = Y + n * dT4Y;
  center_2 = Y + n * (-dT4Y);
  return true;
}


double angle(Coor a, Coor b, Coor n) {
  Coor c;
  double s, p;
  a = norm(a);
  b = norm(b);
  s = a * b;
  c = a % b;
  p = c * n;
  if (p < 0) {
    double r = 2 * M_PI - acos(s);
    return std::isnan(r) ? 1e-02 : r;
  }
  else {
    double r = acos(s);
    return std::isnan(r) ? 1e-02 : r;
  }
}

Coor rotate_vector(Coor x_vec, Coor y_vec, Coor X, double angle, double r) {
  return X + x_vec * (r * cos(angle)) + y_vec * (r * sin(angle));
}

Coor calculate_f2(Atom at, Coor probe_center, double p_r) {
  /*
    Calculate the point f, where the probe sphere and the atom are in contact.
  */
  Coor e;
  e = probe_center - at.crd;
  return at.crd + e * (at.r / (at.r + p_r));
}

double line_dist(Coor t1, Coor e, Coor t2, Coor f) {
  Coor n = norm(e % f);
  Coor t2t1 = t2 - t1;
  return n * t2t1;
}

double line_point_dist(Coor t1, Coor e, Coor t2) {
  Coor X = t1 + project(t2 - t1, e);
  return dist(X, t2);
}
    
bool check_direction(Atom at1, Atom at2, Atom at3, Coor center) {
  /*
    The triangle between the atom centers at1, at2, at3 can be 
    clockwise or anti-clockwise. This is checked here by calculating normal
    to the plane of atoms at1, at2, at3, and a vector from at1.crd to center. If 
    the dot product of these two vectors is greater than zero (the angle between 
    these two vectors must be < pi/2) then we return true else false.
  */
  Coor a21, a31, center_at1;
  Coor n;
  a21 = at2.crd - at1.crd;
  a31 = at3.crd - at1.crd;
  n = a21 % a31;
  center_at1 = center - at1.crd;
  //  cout << "n*center_at1 = " << n * center_at1 << endl; 
  if (n * center_at1 > 0) 
    return true;
  else 
    return false;
}
