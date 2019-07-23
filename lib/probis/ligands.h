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

#ifndef _LIGANDS_H
#define _LIGANDS_H

#include "const.h"
#include "ligand.h"
#include "residue.h"

class Molecule;
class Coor;
class BSite;
class Residue;
class Grid;

class Ligands {
 public:
 Ligands() : nPP(0), nSL(0), nNU(0), nIO(0), ncolumns(3), znPP(0), znSL(0), znNU(0), znIO(0) {}
  ~Ligands();
  void read(Molecule*);
  void generate(Molecule*);
  void center(Molecule*);
  void cluster();
  void output();
  void init_grid_ligands(Grid*, float);
  void zip_points();
  void initialize_data();
  void pdb_header(const string);
  void write_rotate_ligands(Molecule*, bool);
 private:
  bool check_cluster(double, int*, int, int, double**, Coor&);
  void do_ligands(int, int, double, int*, int*, double**, Ligand**, bstype);

  multimap<Residue*, Ligand*, rescomp> ligands;
  set<Ligand*, ligcomp> lig_query; // povezuje ligand (Ligand) z query aminokislino (ChResi)
  vector<BSite*> bsite;    // binding site-i za vse tipe ligandov
  int nPP, nSL, nNU, nIO;
  const int ncolumns;
  double **dataPP, **dataSL, **dataNU, **dataIO;
  Ligand **pligPP, **pligSL, **pligNU, **pligIO;

  int znPP, znSL, znNU, znIO;
  double **zdataPP, **zdataSL, **zdataNU, **zdataIO;
  Ligand **zpligPP, **zpligSL, **zpligNU, **zpligIO;

};

extern Ligands *lig;



#endif // _LIGANDS_H
