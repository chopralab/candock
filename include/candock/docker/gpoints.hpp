#ifndef GPOINTS_H
#define GPOINTS_H
#include "candock/geometry/geometry.hpp"
#include "candock/helper/array1d.hpp"
#include "candock/helper/array3d.hpp"
#include "candock/helper/error.hpp"
#include "candock/centro/centroids.hpp"

namespace candock {

namespace score {
        class Score;
}

namespace docker {

        class Gpoints {
        public:

                struct IJK {
                        int i, j, k;
                        IJK operator+ (const IJK &right) const {
                                return IJK {i + right.i, j + right.j, k + right.k};
                        }
                        IJK operator- (const IJK &right) const {
                                return IJK {i - right.i, j - right.j, k - right.k};
                        }
                        friend ostream &operator<< (ostream &os, const IJK &ijk) {
                                os << ijk.i << " " << ijk.j << " " << ijk.k << endl;
                                return os;
                        }
                };

                struct Gpoint {
                        geometry::Point __crd;
                        IJK __ijk;
                        Array1d<double> __energy; // ligand idatm type to energy
                        geometry::Point &crd() {
                                return __crd;
                        }
                        double energy (const int l) const {
                                return __energy.data[l];
                        }
                        const geometry::Point &crd() const {
                                return __crd;
                        }
                        IJK &ijk() {
                                return __ijk;
                        }
                        const IJK &ijk() const {
                                return __ijk;
                        }
                        void distance (double) const {} // just dummy : needed by grid
                };

                typedef vector<Gpoint> GpointVec;
                typedef vector<Gpoint *> PGpointVec;

        private:

                map<int, GpointVec> __gridpoints;
                map<int, Array3d<Gpoint *>> __gmap;
                const score::Score *__score;
                const set<int> *__ligand_idatm_types;


                void __identify_gridpoints (const double &grid_spacing, const double &radial_check);
                void __identify_gridpoints (const centro::Centroids &centroids, const molib::Atom::Grid &grid,
                                            const double &grid_spacing, const int &dist_cutoff, const double &excluded_radius,
                                            const double &max_interatomic_distance);

        public:
                Gpoints (const double &grid_spacing, const double &radial_check);
                Gpoints (const score::Score &score, const set<int> &ligand_idatm_types,
                         const centro::Centroids &centroids, const molib::Atom::Grid &grid,
                         const double &grid_spacing, const int &dist_cutoff,
                         const double &excluded_radius, const double &max_interatomic_distance);
                Gpoints (const centro::Centroids &centroids, molib::Atom::Grid &grid,
                         const double &grid_spacing, const int &dist_cutoff,
                         const double &excluded_radius, const double &max_interatomic_distance);
                GpointVec &get_gridpoints0() {
                        if (__gridpoints.count (0) != 0) {
                                return __gridpoints.at (0);
                        }

                        throw Error ("die : no gridpoints0 ?");
                }
                const GpointVec &get_gridpoints0() const {
                        if (__gridpoints.count (0) != 0) {
                                return __gridpoints.at (0);
                        }

                        throw Error ("die : no gridpoints0 ?");
                }
                const map<int, GpointVec> &get_gridpoints() const {
                        return __gridpoints;
                }
                const Gpoint &get_center_point() const;
                const Array3d<Gpoint *> &get_gmap (const int bsite_id) const {
                        if (__gmap.count (bsite_id) != 0) {
                                return __gmap.at (bsite_id);
                        }

                        throw Error ("die : cannot get gmap for bsite #" + std::to_string (bsite_id));
                }

                friend ostream &operator<< (ostream &os, const Gpoints &gpoints);
        };

        ostream &operator<< (ostream &os, const docker::Gpoints::GpointVec &points);
        ostream &operator<< (ostream &os, const docker::Gpoints::PGpointVec &points);
};

}

#endif
