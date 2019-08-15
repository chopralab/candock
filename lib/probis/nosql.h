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
