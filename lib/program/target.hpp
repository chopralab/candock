#ifndef TARGET_H
#define TARGET_H

#include "fragmentligands.hpp"
#include "design/design.hpp"

#include "docker/gpoints.hpp"
#include "molib/molecules.hpp"

#include "modeler/forcefield.hpp"
#include "program/findcentroids.hpp"
#include "program/dockfragments.hpp"
#include "program/linkfragments.hpp"

#include <string>

namespace Program {

        class CANDOCK_EXPORT Target {

                std::unique_ptr <Molib::Molecule>         __protein;
                std::unique_ptr <Score::Score>            __score;
                std::unique_ptr <OMMIface::ForceField>    __ffield;
                std::unique_ptr <Molib::Atom::Grid>       __gridrec;
                std::unique_ptr <Program::FindCentroids>  __centroids;
                std::unique_ptr <Program::DockFragments>  __prepseeds;
                std::unique_ptr <Program::LinkFragments>  __dockedlig;

                void __initialize_score(const FragmentLigands &ligand_fragments);
                void __initialize_ffield();
                void __initialize_kbforce();

        public:
                Target (const std::string &input_name);
				//Target(const Target &) {};
				//Target &operator= (const Target &) {};

                std::set<int> get_idatm_types (const std::set<int> &previous = std::set<int>()) const;

                void find_centroids ();
                void make_gridhcp   (const FragmentLigands &ligand_fragments);
                void dock_fragments (const FragmentLigands &ligand_fragments);
                void link_fragments (const FragmentLigands &ligand_fragments);
                void link_fragments (const Molib::Molecules &ligand_fragments);
                void make_scaffolds (const std::set<std::string>& seeds_to_add, Molib::Molecules& all_designs_out );
                void design_ligands (const std::set<std::string>& seeds_to_add, Molib::Molecules& all_designs_out );

                const DockFragments& get_seeds() const {
                        return *__prepseeds;
                }

                const std::string& name() const {
                        return __protein->name();
                }
        };

}

#endif // TARGET_H
