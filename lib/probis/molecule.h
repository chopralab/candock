#ifndef _MOLECULE_H
#define _MOLECULE_H

#include "const.h"
#include "geo.h"
#include "residue.h"
#include "ligand.h"

class Atom;
class Probe;
class Descriptor;
class Item;
class Subgraph;
class Residue;
class Ligand;
class BSite;
class Grid;
class BioUnit;
class EElement;
class ClusterData;
class Motif;




/* napovemo primerjalne funkcije */
bool by_bstype_bsscore(BSite*, BSite*);
bool by_flx_cluster_score(const pair<Residue*, Residue*>&, const pair<Residue*, Residue*>&);
bool by_noflx_cluster_score(const pair<Residue*, Residue*>&, const pair<Residue*, Residue*>&);
bool by_conservation(const pair<ChResi,int>&, const pair<ChResi, int>&);

class Molecule {
 public:
  Molecule (string, string="", int=0);
  ~Molecule ();
  string name;
  string pdb_id;
  int main_id;  // vsako prilegano strukturo v lista_molekul lahko identificiramo po main_id
  string chain_id;
  int unique_seq_size; // prilegana dolzina (ker je seq tipa multimap se klastri lahko prekrivajo)
  float first_cluster_score; // za prilegani protein, cluster_score prvega klastra

  multimap<Residue*, Residue*, rescomp> align;
  set<Residue*, rescomp> sequence;
  multimap<pair<char,int>, char> imap;
  vector<Molecule*> lista_molekul;   // ujemajoci se deli zaporedij proteinov, s katerimi smo primerjali nas protein
  vector<string> plain_pdb; // pdb oziroma pdbqt file v tekstovni obliki

  map<pair<int, int>, ClusterData* > align_score; // OBSOLETE (ONLY FOR MYSQL) za prilegane proteine mapira (main_id, clus_id) --> ClusterData 
  map<pair<int,int>, BioUnit*> biounit;  //bu je dolocen s stevilko modela in stevilko bu-ja

  Atom *atom, *comp;
  int size;
//  int read_PDB(bool, int, bool, int=_noSaveMem);
  void read_PDB(bool, int, bool, int=_noSaveMem);
  void read_plain_pdb(bool);
  void fn_important_residues();
  void output_rotated_asym(Molecule*, Item*, int);
  void output_rotate_pdb(Molecule*, Item*);
  void output_cons_pdb();
//  void output_json_cons();
  void output_iface_pdb();
  void interface(float);
  void interface_models(float);
  void output_interface_models();
  void append_coor(string="");
  void append_coor_delete_comp(string="");
  void init_grid_molecule(Grid*, float);
  void mark_motif_atoms(Motif*);
  void mark_binding_site(Molecule*, string);
  void all_triples(Probe*&);
  void enumerate_clefts(Probe*);
#ifdef CILE
  void patches_cile(Grid*, Probe*);
#endif
  void output_surface();
  void surface_atoms(Grid*, Probe*);
  void set_descriptors(Descriptor*&, bool);
//  void output(bool=false);
  void output();
  void restore_atoms();
  void restore_descriptors(Descriptor*&, bool);
  void restore_probes(Probe*&);
  void get_seq();
  void binarize_cons_seq();
//  void out_seq();
//  void out_ali_seq();
//  void out_html(Subgraph*);
  void out_info();
//  void out_cons();
  void out_query();
  void out_alignments(Subgraph*, const string&);
  void goodness_of_prediction();
  void check_same_coor(Atom*);
  void fingerprint();
//  void signature_residues2();
//  void conservation0(float);
//  void conservation2(float);
  void conservation(float);
  void statistics();
  pair<Coor, Coor> minmax_coor();
  void set_all_atoms_surface();
  void gen_cons_biounit();
  void remark_pdb(Item*);
  int get_lcolor() { return lcolor; }
 private:
  int nfp;
  int lcolor;
  int surf_size_color[1000];
  void expand(Probe*, int);
#ifdef CILE
  void expand_cile(Probe*, Probe*);
  void dijkstra(set<Probe*>, Probe*);
  void surface_atoms_cile(Grid*, Probe*);
#endif
  void moves_from_centers(Atom*, Atom*, EElement*);
  bool check_next_resi(size_t);
  string remark_scores(Item*);
           
  map<int, string> modelName;
  set<Ligand, ligcomp2> modres;
  map<Ligand, string, ligcomp2> site;
  multimap<Ligand, Ligand, ligcomp2> imap_models;  // kadar gledamo binding site na query-ju do vseh prileganih ligandov (vsak je v svojem MODEL)

};

extern Molecule *mol1, *mol2;


#endif // _MOLECULE_H
