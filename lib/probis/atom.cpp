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

#include "atom.h"
#include <iostream>
#include "eelement.h"
#include "geo.h"

Atom& Atom::operator=(const Atom& other) {
#ifdef VERB
    cout << "Atom:: & operator = (const Atom &) " << endl;
#endif
    if (this != &other) {
        strcpy(tag, other.tag);
        strcpy(resn, other.resn);
        chain_id = other.chain_id;
        resi = other.resi;
        acid_start = other.acid_start;
        crd = other.crd;
        r = other.r;
        het = other.het;
        model = other.model;
    }
    return *this;
}

#ifdef CILE
void Atom::pdb() {
    char line[SMALL];
    sprintf(line, "HETATM%5d%4s%5s%2s%4d%12.3f%8.3f%8.3f%6s%6.2f", this->num,
            this->tag, this->resn, to_string(this->chain_id).c_str(),
            this->resi, this->crd.x, this->crd.y, this->crd.z, "1.00", 1.0);
    cout << line << endl;
}
#endif

void Atom::print_atom(ostream& os) {
    char* buffer = new char[100];
    os << "A>";
    print_sphere(os);
    //  if (!strcmp(chain_id, " ")) sprintf(chain_id, "%s", "0");
    //  sprintf(buffer, "%3d%1s%-4s%4s%5d%2c%3d", helix, " ", tag, resn, resi,
    //  chain_id, colors.size);
    sprintf(buffer, "%3d%1s%-4s%4s%5d%2c%3d", model, " ", tag, resn, resi,
            chain_id, colors.size);
    os << buffer;
    for (int i = 0; i < colors.size; i++) {
        sprintf(buffer, "%5d", colors.color[i]);
        os << buffer;
    }
    os << endl;
    delete[] buffer;
}

bool Atom::overlap(Coor center) {
    for (EElement* tmp = neighb; tmp != NULL; tmp = (EElement*)tmp->next)
        if (tmp->item->r + PROBE - dist(tmp->item->crd, center) > EPS)
            return true;
    return false;
}

Atom::~Atom() { delete_neighbor_list(); }

void Atom::delete_neighbor_list() {
    EElement* tmp;
    while (neighb != NULL) {
        tmp = neighb;
        neighb = (EElement*)neighb->next;
        delete tmp;
    }
    neighb = NULL;
}

void Atom::check_aasize() {
/*
  Check if aminoacid is of right size. If not, write -1 into aasize. This is
  done only for the first atom in each aa,
  e.g., the one to which acid_start points.
  Note: resn and resi must be determined before calling check_aasize.
*/
#ifdef VERB
    cout << "Atom::check_aasize() aasize = " << aasize << endl;
#endif

    if (!strncmp(resn, "ALA", 3)) {
        if (aasize < 5) aasize = -1;
    } else if (!strncmp(resn, "ARG", 3)) {
        if (aasize < 11) aasize = -1;
    } else if (!strncmp(resn, "ASN", 3)) {
        if (aasize < 8) aasize = -1;
    } else if (!strncmp(resn, "ASP", 3)) {
        if (aasize < 8) aasize = -1;
    } else if (!strncmp(resn, "CYS", 3)) {
        if (aasize < 6) aasize = -1;
    } else if (!strncmp(resn, "GLN", 3)) {
        if (aasize < 9) aasize = -1;
    } else if (!strncmp(resn, "GLU", 3)) {
        if (aasize < 9) aasize = -1;
    } else if (!strncmp(resn, "GLY", 3)) {
        if (aasize < 4) aasize = -1;
    } else if (!strncmp(resn, "HIS", 3)) {
        if (aasize < 10) aasize = -1;
    } else if (!strncmp(resn, "ILE", 3)) {
        if (aasize < 8) aasize = -1;
    } else if (!strncmp(resn, "LEU", 3)) {
        if (aasize < 8) aasize = -1;
    } else if (!strncmp(resn, "LYS", 3)) {
        if (aasize < 9) aasize = -1;
    } else if (!strncmp(resn, "MET", 3)) {
        if (aasize < 8) aasize = -1;
    } else if (!strncmp(resn, "PHE", 3)) {
        if (aasize < 11) aasize = -1;
    } else if (!strncmp(resn, "PRO", 3)) {
        if (aasize < 7) aasize = -1;
    } else if (!strncmp(resn, "SER", 3)) {
        if (aasize < 6) aasize = -1;
    } else if (!strncmp(resn, "THR", 3)) {
        if (aasize < 7) aasize = -1;
    } else if (!strncmp(resn, "TRP", 3)) {
        if (aasize < 14) aasize = -1;
    } else if (!strncmp(resn, "TYR", 3)) {
        if (aasize < 12) aasize = -1;
    } else if (!strncmp(resn, "VAL", 3)) {
        if (aasize < 7) aasize = -1;
    }
    if (aasize == -1)
        cout << "Warning (ATOM) : Problems with residue " << resn << " " << resi
             << endl;
}
