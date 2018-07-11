#include "probe.h"
#include "debug.hpp"

void Probe::output(int lcolor, bool _motif) {
    /*
      Izpisemo probe centre. Ce je _motif == true, potem izpisemo samo tiste, ki
      pripadajo najvecjemu surface-u.
      Ce na primer izberemo kot motif dva patcha na povrsini, ki nista skupaj,
      potem bi v nasprotnem primeru imeli
      dva surface-a z dvema razlicnima color.
     */
    for (Probe* p = this; p != NULL; p = p->next) {
        if (_motif && lcolor != p->color) continue;
        p->print_probe();
    }
}

void Probe::invert(Probe*& first) {
    Probe* d = first;
    Probe* prev = NULL;
    while (d) {
        Probe* next = d->next;
        d->next = prev;
        prev = d;
        d = next;
    }
    first = prev;
}

void Probe::print_probe(ostream& os) {
    char* buffer = new char[100];
    os << "P>";
    print_sphere(os);
    sprintf(buffer, "%5d", color);
    os << buffer << endl;
    delete[] buffer;
}

#ifdef CILE
bool by_probe_dist(const Probe* i, const Probe* j) { return i->dist < j->dist; }

// void Probe::delete_neighbor_list() {
//  EElement *tmp;
//  while(neighb!=NULL) {
//    tmp = neighb;
//    neighb = (EElement*) neighb->next;
//    delete tmp;
//  }
//  neighb = NULL;
//}
//
// void Probe::init_grid(Grid *grid, Molecule *m, float distance) {
//  /*
//     V probe.neighb zapisemo bliznje atome.
//  */
//
//  for (Probe *p = this; p != NULL; p=p->next) {
//    p->delete_neighbor_list();
//  }
//
//  /* v gridu vsi atomi (tudi vsi modeli & hetero), ki so v mejah (0,0,0)
//  (NUM_CELLS,NUM_CELLS,NUM_CELLS) - make_grid sam preveri!!! */
//  for (Atom *a = m->atom; a != NULL; a=a->next) {
//    grid->make_grid(a);
//  }
//
//  /* za vsak probe center naredimo listo sosednjih atomov */
//  for (Probe *p = this; p != NULL; p=p->next) {
//    grid->volume_slice_probe(p, distance + PROBE + MAXR, distance);
//  }
//
//}

void Probe::pdb() {
    char line[SMALL];
    sprintf(line, "HETATM%5d%4s%5s%2s%4d%12.3f%8.3f%8.3f%6s%6.2f", 1, "CX",
            "XXX", "X", this->num, this->crd.x, this->crd.y, this->crd.z,
            "1.00", this->dist);
    cout << line << endl;
}

#endif  // CILE

void Probe::output_probe() {
    /*
      Output probe centers and color them according to which surface they belong
      to.
    */

    char line[300];
    char tmp[20];
    int i = 0;
    Probe* p = this;
    while (p != NULL) {
        sprintf(line,
                "HETATM                                                        "
                "                                  ");
        sprintf(tmp, "%d", i + 1);
        sprintf(&line[11 - strlen(tmp)],
                "%s                                                            "
                "               ",
                tmp);
        sprintf(&line[13], "%s",
                "S   SUR     1                                                 "
                "                       ");
        sprintf(tmp, "%2.3f", p->crd.x);
        sprintf(&line[38 - strlen(tmp)],
                "%s                                                            "
                "               ",
                tmp);
        sprintf(tmp, "%2.3f", p->crd.y);
        sprintf(&line[46 - strlen(tmp)],
                "%s                                                            "
                "               ",
                tmp);
        sprintf(tmp, "%2.3f", p->crd.z);
        sprintf(&line[54 - strlen(tmp)],
                "%s  1.00                                                      "
                "               ",
                tmp);
        sprintf(tmp, "%2.2f", (double)p->color);
        sprintf(&line[66 - strlen(tmp)],
                "%s                                                            "
                "               ",
                tmp);
        line[66] = '\0';
        cout << line << endl;
        p = p->next;
        i++;
    }
    cout << "END" << endl;
}

void Probe::reset_visited() {
    Probe* p = this;
    while (p != NULL) {
        p->visited = false;
        p = p->next;
    }
}

void Probe::free() {
    dbgmsg("Deleting probe atoms ... ");
    Probe *p = this, *tmp;
    int i = 0;
    while (p != NULL) {
        tmp = p;
        p = p->next;
        delete tmp;
        i++;
    }
    dbgmsg("..deleted " << i << " probe atoms.");
}
