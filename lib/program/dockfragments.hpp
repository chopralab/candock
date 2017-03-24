#ifndef DOCKFRAGMENTS_H
#define DOCKFRAGMENTS_H

#include "programstep.hpp"

#include "findcentroids.hpp"
#include "fragmentligands.hpp"

#include "score/score.hpp"
#include "docker/dock.hpp"
#include "pdbreader/nrset.hpp"

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

                std::vector<std::pair<double, std::string>> get_best_seeds () const;

        };

}

#endif // DOCKFRAGMENTSSTEP_H
