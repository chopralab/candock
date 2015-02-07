#ifndef _LIGAND_H
#define _LIGAND_H

#include "const.h"
#include "sphere.h"

class ClusterData;
class Element;
class Ligand;

/* napovemo primerjalne funkcije */
bool by_type_size_neighb(Ligand*, Ligand*);
bool by_cluster_score(Ligand*, Ligand*);
//bool by_molecule_clid_model_chid_resi(Ligand*, Ligand*);
bool by_ligcomp(Ligand*, Ligand*); // NOVO

class Ligand : public Sphere {
 public:
 Ligand() : Sphere(0), model(0), pdb_id(""), acid(0), chain_id(' '), resn(""), resi(0), cd(NULL), type(_pp), size_neighb(0), neighb(NULL) {}
// Ligand(int M, string P, char CH, string LC, int RI, ClusterData *CD) : Sphere(0), model(M), pdb_id(P), chain_id(CH), resn(LC), resi(RI), cd(CD), type(_pp), size_neighb(0), neighb(NULL) {}
 Ligand(int M, string P, char ACH, char CH, string LC, int RI, ClusterData *CD) : Sphere(0), model(M), pdb_id(P), acid(ACH), chain_id(CH), resn(LC), resi(RI), cd(CD), type(_pp), size_neighb(0), neighb(NULL) {}  // NOVO
 Ligand(int M, string P, char CH, string LC, int RI) : Sphere(0), model(M), pdb_id(P), acid(0), chain_id(CH), resn(LC), resi(RI), cd(NULL), type(_pp), size_neighb(0), neighb(NULL) {}
// Ligand(int M, string P, char CH, string LC, int RI, bstype T) : Sphere(0), model(M), pdb_id(P), chain_id(CH), resn(LC), resi(RI), cd(NULL), type(T), size_neighb(0), neighb(NULL) {}
 Ligand(int M, string P, char ACH, char CH, string LC, int RI, bstype T) : Sphere(0), model(M), pdb_id(P), acid(ACH), chain_id(CH), resn(LC), resi(RI), cd(NULL), type(T), size_neighb(0), neighb(NULL) {}  // NOVO
 Ligand(int M, char CH, int RI) : Sphere(0), model(M), pdb_id(""), acid(0), chain_id(CH), resn(""), resi(RI), cd(NULL), type(_pp), size_neighb(0), neighb(NULL) {}
  ~Ligand();
  int model;
  string pdb_id;
  char acid;  // NOVO
  char chain_id;
  string resn;
  int resi;   // residue number za ligand (resn) v PDB datoteki
  ClusterData *cd;
  set<ChResi> rlist; // set bliznjih query residuejev (vsi pripadajo istemu cluster_id)
  bstype type;
  int size_neighb;
  Element *neighb;
  Ligand* first_neighb();
  Ligand* next_neighb();
  void delete_neighb();
  bool operator<(const Ligand &other) const;
//#ifdef VERB
//  void output() { cout << "model = " << model << " pdb_id = " << pdb_id << " chain_id = " << chain_id << " resn = " << resn << " resi = " << resi << " cd = " << cd << endl; }
  void output() { cout << "model = " << model << " pdb_id = " << pdb_id << " acid = " << acid << " chain_id = " << chain_id << " resn = " << resn << " resi = " << resi << " cd = " << cd << endl; }  // NOVO
//#endif
 private:
  void free_neighb();
  Element *curr, *prev;
};

/* za primerjavo Ligands (cd ne sme biti NULL!)
resi   pdb     ch model    resn  cluster_id
---------------------------------------------
197    1allA    A     1    ATP       0
197    1allA    A     1    ATP       1
197    1allA    A     1    HOH       1
197    1j3gA    A     2    ATP       0
*/
struct ligcomp {
  bool operator() (const Ligand*, const Ligand*) const;
};

/* za primerjavo Ligands (ce gre za lvalue in ni dolocen cd->cluster_id!)
resi   pdb     ch model    resn
---------------------------------
197    1allA    A     1    ATP   
197    1allA    A     1    HOH   
197    1j3gA    A     2    ATP   
*/
struct ligcomp2 {
  bool operator() (const Ligand&, const Ligand&) const;
};


struct ligcomp_bsite {
  bool operator() (const Ligand*, const Ligand*) const;
};






#endif // _LIGAND_H
