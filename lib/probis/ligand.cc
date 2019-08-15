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

#include "ligand.h"
#include "molecule.h"
#include "eelement.h"
#include "clusterdata.h"

bool by_type_size_neighb(Ligand *i, Ligand *j) { 
  if (i->type < j->type)
    return true;
  else if (i->type == j->type)
    return i->size_neighb > j->size_neighb;
  return false;
    }

bool by_cluster_score(Ligand *i, Ligand *j) { return i->cd->cluster_score < j->cd->cluster_score; }



Ligand::~Ligand() {
#ifdef VERB
  cout << "Ligand::~Ligand  :  zbrisem neighb!!!" << endl;
#endif
  free_neighb();
}

void Ligand::free_neighb() {
  Element *tmp;
  while(neighb!=NULL) {
    tmp = neighb;
    neighb = neighb->next;
    delete tmp;
  }
  neighb = NULL;
}

Ligand* Ligand::first_neighb() {
  curr = neighb;
  prev = NULL;
  return curr==NULL ? NULL : (Ligand*)curr->item;
}

Ligand* Ligand::next_neighb() {
  if (curr) {
    prev = curr;
    curr = curr->next;
  }
  return curr==NULL ? NULL : (Ligand*)curr->item;
}

void Ligand::delete_neighb() {
  if (prev) {
    prev->next = curr->next;
    delete curr;
    curr = prev->next;
  }
  else {
    neighb = curr->next;
    delete curr;
    curr = neighb;
  }
}

inline bool Ligand::operator<(const Ligand &other) const {
  if (this->cd < other.cd) {
    return true;
  }
  else if (this->cd == other.cd) {
    if (this->model < other.model) {
      return true;
    }
    else if (this->model == other.model) {
      if (this->chain_id < other.chain_id) {
        return true;
      }
      else if (this->chain_id == other.chain_id) {
        if (this->resi < other.resi) {
          return true;
        }
        else if (this->resi == other.resi) {
          if (this->pdb_id.compare(other.pdb_id) < 0) {
            return true;
          }
          else if (this->pdb_id.compare(other.pdb_id) == 0) {
            if (this->acid < other.acid) {  // NOVO
              return true;
            }            
            else if (this->acid == other.acid) {
              if (this->resn.compare(other.resn) < 0) {
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}


bool ligcomp::operator() (const Ligand *lhs, const Ligand *rhs) const
{
  return *lhs < *rhs;
}


bool ligcomp2::operator() (const Ligand &lhs, const Ligand &rhs) const
{
  return lhs < rhs;
}


bool ligcomp_bsite::operator() (const Ligand *lhs, const Ligand *rhs) const
{
  // najprej sortiramo po cluster_score in stevilu binding site residue-jev
  if (lhs->cd->cluster_score > rhs->cd->cluster_score) {
    return true;
  }
  else if (lhs->cd->cluster_score == rhs->cd->cluster_score) {
    if (lhs->rlist.size() > rhs->rlist.size()) {
      return true;
    }
    else if (lhs->rlist.size() == rhs->rlist.size()) {
      return *lhs < *rhs;
    }
  }
  return false;
}



bool by_ligcomp(Ligand *lhs, Ligand *rhs) {  // NOVO
  return *lhs < *rhs;
}

