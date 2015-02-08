#include "desc.h"
#include "atom.h"
#include "eelement.h"
#include "grid.h"
#include "probe.h"
#include "product.h"
#include "debug.hpp"

void Descriptor::mark_backbone() {
  /* 
     Oznacimo deskriptorje (desc->bb = true), ki pripadajo backbone aminokislinam. 
  */

  Descriptor *d = this;

  while (d) {
    /* ce je backbone, deskriptor oznacimo */
    if(strcmp(d->atom->tag, "N") == 0 || strcmp(d->atom->tag, "O") == 0) { //PEP 
//      dbgmsg("BB>"<< d->atom->tag );
      d->bb = true;
    }
    else 
      d->bb = false;
    
    d = d->next;
  }
}

void Descriptor::reset_visited() {
 Descriptor *d = this; 
 while(d != NULL) { 
   d->visited = false; 
   d = d->next; 
 }
}

void Descriptor::output() {
  Descriptor *d = this;
  while (d) {
    d->print_desc();
    d = d->next;
  }
}

void Descriptor::invert(Descriptor*& first) {
  Descriptor *d = first;
  Descriptor *prev = NULL;
  while (d) {
    Descriptor *next = d->next;
    d->next = prev;
    prev = d;
    d = next;
  }
  first = prev;
}


void Descriptor::print_desc(ostream &os) {
  char *buffer = new char[100];
  os << "D>";
  print_sphere(os);
  sprintf(buffer, "%3d%3d%5d%3c", s->mnsp, psurf, atom->num, atom->chain_id);
  os << buffer << endl;
  delete buffer;
}

//void Descriptor::trim(string chain1, string chain2) {
//  /*
//    We take only descriptor of chain1 which is closer than PCUT 
//    to some atom in chain2.
//  */
//  Descriptor *d = this;
//  while (d != NULL) {
//    d->psurf = 0;
//    if (chain1.find(d->atom->chain_id) != string::npos) 
//      for (EElement *tmp = d->atom->neighb; tmp != NULL; tmp = (EElement*) tmp->next) {
//        if ( chain2.find(((Atom*)tmp->item)->chain_id) != string::npos )
//          if (dist(tmp->item->crd, d->atom->crd) < PCUT) {
//            d->psurf = 1;
//          }
//      }
//    d = d->next;
//  }
//}

Descriptor::~Descriptor() {
  delete s;  
  free_neighb();
}

void Descriptor::calculate_sum() {
  /*
    The sum of the vectors from descriptors to the closest probe centers (in the neighb array) is calculated.
    This sum vector is scaled with the number of probe centers.
  */
  Descriptor *d = this;
  while (d != NULL) {
    set_zero(d->sum);
    for (Element *tmp=d->neighb; tmp!=NULL; tmp=tmp->next)
      d->sum = d->sum + tmp->item->crd - d->crd;
    d->sum = norm(d->sum);
    d->sum = d->crd + d->sum;
    d = d->next;
  }
}

//    void Descriptor::set_schmitt_weights() {
//      /*
//        Each schmitt descriptor is assigned its neighbours less than SCHMITT distance away. 
//        The result is in desc.s.C[j][i] array, j indicating the descriptor tag (0 = al, 1 = pi, 2 = ac, 
//        3 = do, acdo are assigned to 2 and 3), and i is the indicator of distance from the central descriptor 
//        (see defined variable F = 4 means 0.25 Angstrom precision. Second part smoothes the distribution (0,1,0) 
//        becomes (0.125,0.75,0.125), etc.
//    
//        DEC/06/2009 OBSOLETE!
//      */
//      Descriptor *d = this;
//      double f1 = 0.75, f2 = 0.125;
//      double ci = 0.0, cj = 0.0, cj_prev = 0.0;
//      double distance;
//      int i, j;
//      while (d != NULL) {
//        for (Element *tmp=d->neighb; tmp!=NULL; tmp=tmp->next) {
//          if ((distance = dist(d->crd, tmp->item->crd)) < SCHMITT) {
//            switch (((Descriptor*) tmp->item)->s->mnsp) 
//              { 
//              case(0x01): d->s->C[0][(int) (distance*F) + 2]++; break;
//              case(0x02): d->s->C[1][(int) (distance*F) + 2]++; break;
//              case(0x04): d->s->C[2][(int) (distance*F) + 2]++; break;
//              case(0x08): d->s->C[3][(int) (distance*F) + 2]++; break;
//              case(0x0C): d->s->C[2][(int) (distance*F) + 2]++;d->s->C[3][(int) (distance*F) + 2]++; break;
//                //          case(0x0C): d->s->C[4][(int) (distance*F) + 2]++; break;
//              }
//          }
//        }
//        //    output_descriptor(d, NULL, 0);
//        for (j = 0; j < SCHMITT_NUM; j++) {
//          cj_prev = 0.0;
//          for (i = 2; i <= (int) ceil(SCHMITT*F + 1); i++) {
//            ci = f1*d->s->C[j][i];
//            cj = f2*d->s->C[j][i];
//            d->s->C[j][i-1]+=cj;
//            d->s->C[j][i] = ci + cj_prev;
//            cj_prev = cj;
//          }
//          d->s->C[j][i] = cj;
//        }
//        //    for (int i=0; i < SCHMITT_NUM; i++) { for (int j=0; j < SCHMITT*F; j++) dbgmsg(d->s->C[i][j] << " "; cout ); }
//        //    output_descriptor(d, NULL, 0);
//        d = d->next;
//      }
//    }

void Descriptor::init_grid_descriptors(Grid *grid, double max_dist, double distance) {
  /* 
     Descriptor is a child class of Atom, therefore all Atom-s methods and variables can be used.
     Here the descriptor coordinates are mapped to box-shaped segments. First, the descriptors with 
     the minimum x-coordinate is calculated. Second, the descriptors are mapped to grid according
     to their relative coordinates.
     NOTE: OCT/22/2008 init_grid_probe can be called before this (state4, database search), therefore
     neighb list has to be deleted for each descriptor!

     Input is desc.crd;
     Output is grid;
     
     Next, the list of neighboring descriptors, less than CUTOFF_SECOND Angstroms apart, based on grid is generated.
     Input is grid;
     Output is desc.neighb;
  */
  Descriptor *d;

  d = this;
  while (d != NULL) {
    grid->make_grid(d);
    d = d->next;
  }
  d = this;
  while (d != NULL) {
    d->free_neighb();
    grid->volume_slice_desc(d, max_dist, distance);
    d = d->next;
  }
}

void Descriptor::init_grid_probe(Grid *grid, int lcolor, Probe *probe, double max_dist, double distance) {
  /* 
     NOTE: before this procedure, init_grid_descriptors was called which filled the neighb list. 
     Therefore the neighb list has to be deleted for each descriptor.
     Probe is a child class of Atom, therefore all Atom-s methods and variables can be used.
     Here the probe coordinates are mapped to box-shaped segments. First, the probe with 
     the minimum x-coordinate is calculated. Second, the probes are mapped to grid according
     to their relative coordinates.
     Input is desc.crd;
     Output is grid;
     
     For each descriptor, the list of neighboring probes, less than SURF Angstroms apart, 
     based on grid is generated.
     Input is grid, desc;
     Output is desc->neighb;
  */
  Probe *p = probe;
  Descriptor *d = this;
//  grid->init_min_crd(p->crd);
//  grid->init_min_crd();
//  while (p != NULL) {
//    grid->set_min_crd(p);
//    p = p->next;
//  }
  p = probe;
  while (p != NULL) {
    // uplostevamo samo tiste probe, ki so na najvecji povrsini
    if (p->color == lcolor) {
      grid->make_grid(p);
    }
    p = p->next;
  }
  while (d != NULL) {
    d->free_neighb();
    grid->volume_slice_desc(d, max_dist, distance); 
    d = d->next;
  }
}

void Descriptor::free_neighb() {
  Element *tmp;
  while(neighb!=NULL) {
    tmp = neighb;
    neighb = neighb->next;
    delete tmp;
  }
  neighb = NULL;
}

void Descriptor::free() {
  dbgmsg("Deleting descriptors ..." );
  Descriptor *d, *desc;
  Column *c;
  int i = 0;
  desc = this;
  while (desc != NULL) {
    d = desc;
    while (desc->column != NULL) {
      c = desc->column;
      desc->column = desc->column->next;
      delete c;
    }
    desc = desc->next;
    delete d;
    i++;
  }
  dbgmsg("Deleted " << i << " descriptors ..." );
}

void Descriptor::free_columns() {
      dbgmsg("Deleting only columns of descriptors (not descriptors)..." );
      Column *c;
      Descriptor *d = this;
      while (d != NULL) {
        while (d->column != NULL) {
          c = d->column;
          d->column = d->column->next;
          delete c;
        }
        d->column = NULL;
        d->current = NULL;
//        delete d->s;
//        d->s = new Schmitt();
        d = d->next;
      }
}

void Descriptor::output_descriptors() {
  /*
    Output all descriptors colored with color.
  */
  char line[300];
  char tmp[20];
  int i = 0;
  Descriptor *d;
  d = this;
  while (d != NULL) {
    sprintf(line, "HETATM                                                                                          ");
    sprintf(tmp, "%d", i + 1);
    sprintf(&line[11 - strlen(tmp)], "%s                                                                           ", tmp);
    sprintf(&line[13], "%s", "S   DSC     1                                                                        ");
    sprintf(tmp, "%2.3f", d->crd.x);
    sprintf(&line[38 - strlen(tmp)], "%s                                                                           ", tmp);
    sprintf(tmp, "%2.3f", d->crd.y);
    sprintf(&line[46 - strlen(tmp)], "%s                                                                           ", tmp);
    sprintf(tmp, "%2.3f", d->crd.z);
    sprintf(&line[54 - strlen(tmp)], "%s  1.00                                                                     ", tmp);
    //    sprintf(tmp, "%2.2f", (double) d->s->mnsp);
    sprintf(tmp, "%2.2f", (double) d->psurf);
    sprintf(&line[66 - strlen(tmp)], "%s                                                                           ", tmp);
    line[66] = '\0';
    cout << line << endl;
    i++;
    d = d->next;
  }
  cout << "END" << endl;
}

void Descriptor::insert(Row *r) {
  Column *tmp;
  tmp = new Column();
  tmp->item = r;
  tmp->next = column;
  column = tmp;
  current = column;
}

Row* Descriptor::get() {
  Column *tmp;
  if (current != NULL) {
    tmp = current;
    current = current->next;
    return tmp->item;
  }
  else {
    current = column;
    return NULL;
  }
}
