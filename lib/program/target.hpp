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

#ifndef TARGET_H
#define TARGET_H

#include "fragmentligands.hpp"
#include "design/design.hpp"

#include "docker/gpoints.hpp"
#include "molib/molecules.hpp"

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
                                filename(file), protein (rec),
                                score (nullptr), ffield (nullptr),
                                gridrec (nullptr), centroids (nullptr),
                                prepseeds (nullptr), dockedlig (nullptr) {
                        }

                        virtual ~DockedReceptor();

                        DockedReceptor (const DockedReceptor &rhs) = default;
                        DockedReceptor &operator= (const DockedReceptor &rhs) = default;

                        std::string          filename;
                        Molib::Molecule      &protein;
                        Molib::Score         *score;
                        OMMIface::ForceField *ffield;
                        Molib::Atom::Grid    *gridrec;
                        FindCentroids        *centroids;
                        DockFragments        *prepseeds;
                        LinkFragments        *dockedlig;
                };

                Molib::Molecules            __receptors;
                std::vector<DockedReceptor> __preprecs;

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
