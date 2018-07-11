#include "schmitt.h"

bool Schmitt::compare_single(Schmitt* s) {
    if (mnsp & s->mnsp) return true;
    return false;
}
