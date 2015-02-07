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
