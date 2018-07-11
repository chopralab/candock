#include "residue.h"

bool rescomp::operator()(const Residue* lhs, const Residue* rhs) const {
    if (lhs->chain_id < rhs->chain_id) {
        return true;
    } else if (lhs->chain_id == rhs->chain_id) {
        if (lhs->resi < rhs->resi) {
            return true;
        } else if (lhs->resi == rhs->resi) {
            if (lhs->m < rhs->m) {
                //      if (lhs->m && rhs->m && lhs->m < rhs->m) {
                //      if (lhs->cd && rhs->cd && lhs->cd->m < rhs->cd->m) {
                return true;
            }
        }
    }
    return false;
}
