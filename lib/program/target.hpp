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

        // TODO: Implement as a templated_map<Molecule,Target,Target> ????? Or as a Molecules?????

        class FindCentroids;
        class DockFragments;
        class LinkFragments;

        class CANDOCK_EXPORT Target {

                // FIXME: There's a better to design this, but this works for *now*
                // TODO:  Consider using ProgramSteps instead of named things?
                struct CANDOCK_EXPORT DockedReceptor {

                        DockedReceptor (Molib::Molecule &rec, const std::string &file) :
                                filename(file), protein (rec){
                        }

                        DockedReceptor (const DockedReceptor &rhs) = delete;
                        DockedReceptor &operator= (const DockedReceptor &rhs) = delete;

                        std::string          filename;
                        Molib::Molecule      &protein;
                        std::unique_ptr <Molib::Score>            score;
                        std::unique_ptr <OMMIface::ForceField>    ffield;
                        std::unique_ptr <Molib::Atom::Grid>       gridrec;
                        std::unique_ptr <Program::FindCentroids>  centroids;
                        std::unique_ptr <Program::DockFragments>  prepseeds;
                        std::unique_ptr <Program::LinkFragments>  dockedlig;
                };

                Molib::Molecules            __receptors;
                std::vector<std::unique_ptr<DockedReceptor>> __preprecs;

                void __initialize_score(const FragmentLigands &ligand_fragments);
                void __initialize_ffield();
                void __initialize_kbforce();

        public:
                Target (const std::string &input_name);
                Target (const Target &) = delete;
                Target &operator= (const Target &) = delete;

                std::set<int> get_idatm_types (const std::set<int> &previous = std::set<int>());

                // TODO: Instead of named function, pass in fully initiallized ProgramSteps????????
                void find_centroids();
                void make_gridhcp  (const FragmentLigands &ligand_fragments);
                void dock_fragments (const FragmentLigands &ligand_fragments);
                void link_fragments (const FragmentLigands &ligand_fragments);
                void make_scaffolds (FragmentLigands &ligand_fragments, const std::set<std::string> &seeds_to_add);
                void design_ligands (FragmentLigands &ligand_fragments, const std::set<std::string> &seeds_to_add);

                static void make_objective();

                // TODO: Ideally this would be done internally.....
                std::multiset<std::string> determine_overlapping_seeds (const int max_seeds, const int number_of_occurances) const;

                static std::set<std::string> determine_non_overlapping_seeds (const Target &targets, const Target &antitargets);
        };

}

#endif // TARGET_H
