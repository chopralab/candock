#ifndef FILEOUT_H
#define FILEOUT_H

#include <cmath>
#include <iosfwd>
#include "candock/molib/molecule.hpp"

namespace candock {

namespace fileout {
void print_complex_pdb(std::ostream& ss, const molib::Molecule& ligand,
                       const molib::Molecule& receptor, const double energy,
                       const double potential = 0.0, const int model = 1,
                       const size_t max_clq_id = 1,
                       const double rmsd = std::nan(""));

void print_mol2(std::ostream& ss, const molib::Molecule& ligand);
}
}

#endif