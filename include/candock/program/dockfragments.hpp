#ifndef DOCKFRAGMENTS_H
#define DOCKFRAGMENTS_H

#include "candock/program/programstep.hpp"

#include "candock/program/findcentroids.hpp"
#include "candock/program/fragmentligands.hpp"

#include "candock/score/score.hpp"
#include "candock/docker/dock.hpp"
#include "candock/molib/nrset.hpp"

namespace candock {

namespace Program {

        class DockFragments : public ProgramStep
        {
                const FindCentroids& __found_centroids;
                const FragmentLigands& __fragmented_ligands;

                const score::Score& __score;
                const molib::Atom::Grid& __gridrec;

                const std::string& __name;
                std::string __top_seeds_location;

                molib::NRset __all_seeds;

                void __dock_fragment(int start, const docker::Gpoints& gpoints, const docker::Gpoints& gpoints0);
        protected:
                virtual bool __can_read_from_files();
                virtual void __read_from_files();
                virtual void __continue_from_prev();

        public:
                DockFragments ( const FindCentroids& found_centroids,
                                const FragmentLigands& fragmented_ligands,
                                const score::Score& score,
                                const molib::Atom::Grid& gridrec,
                                const std::string& name
                              );

                virtual ~DockFragments(){}

                std::vector<std::pair<double, std::string>> get_best_seeds () const;

                molib::NRset get_negative_seeds(const std::set<std::string> &seeds, const double max_value) const;
                molib::NRset get_top_seeds(const std::set<std::string> &seeds, const double top_percent) const;
                molib::NRset get_seeds(const molib::Molecule      &ligand, const double top_percent) const;

                docker::Gpoints get_gridhcp();
        };

}

}

#endif // DOCKFRAGMENTSSTEP_H
