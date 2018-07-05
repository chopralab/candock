#ifndef LINKFRAGMENTS_H
#define LINKFRAGMENTS_H

#include "candock/program/programstep.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/molib/molecules.hpp"
#include "candock/score/score.hpp"
#include "candock/modeler/forcefield.hpp"
#include "candock/parser/fileparser.hpp"
#include "candock/program/dockfragments.hpp"

namespace candock {

namespace Program {

        class LinkFragments : public ProgramStep
        {
                molib::Molecules __all_top_poses;
                const DockFragments& __seeds_database;

                const molib::Molecule& __receptor;
                const score::Score& __score;
                const OMMIface::ForceField& __ffield;
                const molib::Atom::Grid& __gridrec;

                int __ligand_cnt = 0;

                void __link_ligand (molib::Molecule& ligand);
        protected:
                virtual bool __can_read_from_files ();
                virtual void __read_from_files ();
                virtual void __continue_from_prev ();
        public:
                LinkFragments ( const molib::Molecule& receptor,
                                const score::Score& score,
                                const OMMIface::ForceField& ffield,
                                const DockFragments& seeds_database,
                                const molib::Atom::Grid& gridrec ) :
                                         __seeds_database(seeds_database),
                                        __receptor(receptor), __score(score),
                                        __ffield(ffield), __gridrec(gridrec) {};

                virtual ~LinkFragments() {}
                void link_ligands (const molib::Molecules& ligands);
                const molib::Molecules& top_poses() const { return __all_top_poses; }
                void clear_top_poses() { __all_top_poses.clear(); };
        };
}

}

#endif // LINKFRAGMENTSSTEP_H
