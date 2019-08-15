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

#ifndef DOCKFRAGMENTS_H
#define DOCKFRAGMENTS_H

#include "programstep.hpp"

#include "findcentroids.hpp"
#include "fragmentligands.hpp"

#include "score/score.hpp"
#include "docker/dock.hpp"
#include "molib/nrset.hpp"

namespace Program {

        class DockFragments : public ProgramStep
        {
                const FindCentroids& __found_centroids;
                const FragmentLigands& __fragmented_ligands;

                const Molib::Score& __score;
                const Molib::Atom::Grid& __gridrec;

                const std::string& __name;
                std::string __top_seeds_location;

                Molib::NRset __all_seeds;

                void __dock_fragment(int start, const Docker::Gpoints& gpoints, const Docker::Gpoints& gpoints0);
        protected:
                virtual bool __can_read_from_files();
                virtual void __read_from_files();
                virtual void __continue_from_prev();

        public:
                DockFragments ( const FindCentroids& found_centroids,
                                const FragmentLigands& fragmented_ligands,
                                const Molib::Score& score,
                                const Molib::Atom::Grid& gridrec,
                                const std::string& name
                              );

                virtual ~DockFragments(){}

                std::vector<std::pair<double, std::string>> get_best_seeds () const;

                Molib::NRset get_top_seeds(const std::set<std::string> &seeds, const double top_percent) const;
                Molib::NRset get_top_seeds(const Molib::Molecule      &ligand, const double top_percent) const;
        };

}

#endif // DOCKFRAGMENTSSTEP_H
