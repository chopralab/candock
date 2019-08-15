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

#include "fileout.hpp"

#include "molib/bond.hpp"

namespace Fileout {

        void print_complex_pdb (std::ostream &ss,
                                const Molib::Molecule &ligand, const Molib::Molecule &receptor,
                                const double energy, const double potential,
                                const int model, const size_t max_clq_id,
                                const double rmsd) {

                ss << "MODEL    " << model << "\n";

                ss << "REMARK   1 MINIMIZED COMPLEX OF " << ligand.name() << " AND " << receptor.name() << " WITH SCORE OF " << std::setprecision(6) << energy << endl;

                if (! std::isnan (potential))
                        ss << "REMARK   2 POTENTIAL ENERGY OF " << ligand.name() << " IS " << potential << "\n";

                if (max_clq_id != 0xFFFFFF)
                        ss << "REMARK   3 THIS IS CONFIGUATION " << ligand.name() << " NUMBER " << max_clq_id << "\n";

                if (! std::isnan (rmsd))
                        ss << "REMARK   3 DOCKED CONFORMATION HAS AN RMSD OF " << rmsd << " FROM ORIGINAL STRUCTURE" << "\n";

                int reenum = 0;

                for (auto &patom : receptor.get_atoms()) {
                        ss << std::setw (66) << std::left << *patom;
                        reenum = patom->atom_number();
                }

                for (auto &presd : receptor.get_residues()) {
                        if (help::standard_residues.count (presd->resn()) != 0)
                                continue;

                        // Print internal bonds
                        for (auto pbond : get_bonds_in (presd->get_atoms())) {
                                ss << pbond->get_pdb_remark();
                        }

                        // Print external bonds
                        for (auto pbond : get_bonds_in (presd->get_atoms(), false)) {
                                ss << pbond->get_pdb_remark();
                        }
                }

                ss << "TER" << "\n";

                for (auto &patom : ligand.get_atoms()) {
                        patom->set_atom_number (++reenum);
                        ss << std::setw (66) << std::left << *patom;
                }

                for (auto &pbond : get_bonds_in (ligand.get_atoms())) {
                        ss << pbond->get_pdb_remark();
                }

                ss << "ENDMDL" << "\n";
        }
}
