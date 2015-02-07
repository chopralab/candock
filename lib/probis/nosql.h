#ifndef _NOSQL_H
#define _NOSQL_H

#include "const.h"

class Molecule;
class Item;
class Subgraph;
class Score;

/* napovemo primerjalne funkcije */
bool by_cluster_score_unique_seq_length(Molecule*, Molecule*);

class NoSql {
 public:
	struct WriteStruct {
		string suffix;
		string mol2_pdb_id, mol2_chain_id;
		string mol1_pdb_id, mol1_chain_id;

		string getString(int mol1or2) const {return (mol1or2 == 1 ? mol1_pdb_id + "\t*" + mol1_chain_id : mol2_pdb_id + "\t" + mol2_chain_id) + "\t" + suffix;}
	};
 private:

  void read_file (vector<string>&);
  void read_file (Molecule*, string&);
  void append_file (string);
  void set_db_file(Molecule*);
 public:
	void set_db_file(const string& mol_pdb_id, const string& mol_chain_id);
 private:
  Molecule *mol1;
  string db_file;
 public:
  //~ NoSql(Molecule*);
  //~ NoSql(const string& mol_pdb_id, const string& mol_chain_id);
	NoSql(Molecule *mol, const string &nosql_file) : mol1(mol), db_file(nosql_file) {}
	NoSql(const string &nosql_file) : mol1(0), db_file(nosql_file) {}
  void read();
  void get_rota_trans_resi(Molecule*, Item*);
  void write(Molecule*, Subgraph*, Score*);
	void write(const WriteStruct&);
	void writeStruct(Molecule*, Subgraph*, Score*, WriteStruct&);
};

extern NoSql *nosql;

#endif
