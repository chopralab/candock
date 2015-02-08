#include "eelement.h"

void EElement::exchange(int i, int j) {
  Atom *a;
  Probe *p;

  a = atm[i];
  p = probe[i];

  atm[i] = atm[j];
  probe[i] = probe[j];

  atm[j] = a;
  probe[j] = p;
}
