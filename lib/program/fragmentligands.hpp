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

#ifndef FRAGMENTLIGANDS_H
#define FRAGMENTLIGANDS_H

#include <set>
#include <mutex>

#include "programstep.hpp"
#include "molib/molecules.hpp"
#include "parser/fileparser.hpp"

namespace Program {

        class CANDOCK_EXPORT FragmentLigands : public ProgramStep {
                Molib::Molecules __seeds;
                std::set<int> __ligand_idatm_types;
                std::set<int> __added;

                //Not Used until to ability to avoid rereading from disk is enabled
                //Molib::Molecules __ligands;

                std::mutex __prevent_re_read_mtx;
                std::mutex __add_to_typing_mtx;

                void __fragment_ligands (Parser::FileParser &lpdb, const bool write_out_for_linking, const bool no_rotatable_bond);

        protected:
                virtual bool __can_read_from_files();
                virtual void __read_from_files ();
                virtual void __continue_from_prev ();

        public:
                FragmentLigands() { }
                virtual ~FragmentLigands() {}

                void add_seeds_from_molecules (const Molib::Molecules &molecules);

                const Molib::Molecules &seeds() const {
                        return __seeds;
                }

                const std::set<int> &ligand_idatm_types() const {
                        return __ligand_idatm_types;
                }

                const std::set<int> &all_seeds() const {
                        return __added;
                }

        };

}

#endif
