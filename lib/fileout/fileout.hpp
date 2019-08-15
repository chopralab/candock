/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef FILEOUT_H
#define FILEOUT_H

#include <iosfwd>
#include <cmath>
#include "molib/molecule.hpp"

namespace Fileout {

    void print_complex_pdb( std::ostream &ss,
                        const Molib::Molecule &ligand, const Molib::Molecule &receptor,
                        const double energy, const double potential = 0.0,
                        const int model = 1, const size_t max_clq_id = 1,
                        const double rmsd = std::nan(""));

    void print_mol2( std::ostream &ss, const Molib::Molecule &ligand );
};

#endif