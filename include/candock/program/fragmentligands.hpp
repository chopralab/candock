#ifndef FRAGMENTLIGANDS_H
#define FRAGMENTLIGANDS_H

#include <set>
#include <mutex>

#include "candock/program/programstep.hpp"
#include "candock/molib/molecules.hpp"
#include "candock/parser/fileparser.hpp"

namespace candock {

namespace Program {

        class CANDOCK_EXPORT FragmentLigands : public ProgramStep {
                molib::Molecules __seeds;
                std::set<int> __ligand_idatm_types;
                std::set<int> __added;

                //Not Used until to ability to avoid rereading from disk is enabled
                //molib::Molecules __ligands;

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

                void add_seeds_from_molecules (const molib::Molecules &molecules);
                void add_seeds_from_file      (const std::string &filename);
                void write_seeds_to_file      (const std::string &filename);

                const molib::Molecules &seeds() const {
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

}

#endif
