#ifndef FINDCENTROIDS_H
#define FINDCENTROIDS_H

#include "candock/program/programstep.hpp"
#include "candock/centro/centroids.hpp"

namespace Program {

        class FindCentroids : public ProgramStep {
        protected:
                virtual bool __can_read_from_files();
                virtual void __read_from_files();
                virtual void __continue_from_prev();

                const std::string __filename;
                const std::string __chain_ids;
                const std::string __out_dir;

                Centro::Centroids __result;
                std::string __centroid_file;

        public:
                FindCentroids ( const std::string& filename, const std::string& chain_ids, const std::string& out_dir );
                virtual ~FindCentroids() {}

                const Centro::Centroids& centroids() const {
                        return __result;
                }

        };

}

#endif
