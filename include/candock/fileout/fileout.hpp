#ifndef FILEOUT_H
#define FILEOUT_H

#include <iosfwd>
#include <cmath>
#include "candock/molib/molecule.hpp"

namespace Fileout {

    void print_complex_pdb( std::ostream &ss,
                        const Molib::Molecule &ligand, const Molib::Molecule &receptor,
                        const double energy, const double potential = 0.0,
                        const int model = 1, const size_t max_clq_id = 1,
                        const double rmsd = std::nan(""));

    void print_mol2( std::ostream &ss, const Molib::Molecule &ligand );
};

#endif