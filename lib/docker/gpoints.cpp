/* This is gpoints.cpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

#include "candock/docker/gpoints.hpp"
#include <exception>
#include <iostream>
#include "candock/centro/centroids.hpp"
#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/benchmark.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/molecule.hpp"
#include "statchem/score/score.hpp"

using namespace std;

namespace candock {
namespace docker {

ostream& operator<<(ostream& os,
                    const candock::docker::Gpoints::PGpointVec& points) {
    for (auto& ppoint : points) {
        // to add : output of energies
        os << "ATOM      1   U  DIK     1    " << ppoint->crd().pdb() << endl;
    }
    return os;
}

ostream& operator<<(ostream& os,
                    const candock::docker::Gpoints::GpointVec& points) {
    for (auto& point : points) {
        os << "ATOM      1   U  DIK     1    " << point.crd().pdb() << endl;
    }
    return os;
}

ostream& operator<<(ostream& os, const Gpoints& gpoints) {
    for (auto& kv : gpoints.__gridpoints) {
        const int bsite_id = kv.first;
        os << "MODEL" << setw(9) << right << bsite_id << endl;
        os << kv.second;
        os << "ENDMDL" << endl;
    }
    return os;
}

Gpoints::Gpoints(const double& grid_spacing, const double& radial_check)
    : __score(nullptr), __ligand_idatm_types(nullptr) {
    try {
        __identify_gridpoints(grid_spacing, radial_check);
    } catch (...) {
        log_error << "FAILURE: gridless constructor of Gpoints\n";
        throw;
    }
}
Gpoints::Gpoints(const score::Score& score, const set<int>& ligand_idatm_types,
                 const centro::Centroids& centroids,
                 const molib::Atom::Grid& grid, const double& grid_spacing,
                 const int& dist_cutoff, const double& excluded_radius,
                 const double& max_interatomic_distance)
    : __score(&score), __ligand_idatm_types(&ligand_idatm_types) {
    try {
        __identify_gridpoints(centroids, grid, grid_spacing, dist_cutoff,
                              excluded_radius, max_interatomic_distance);
    } catch (...) {
        log_error << "FAILURE: complete constructor of Gpoints\n";
        throw;
    }
}

Gpoints::Gpoints(const centro::Centroids& centroids, molib::Atom::Grid& grid,
                 const double& grid_spacing, const int& dist_cutoff,
                 const double& excluded_radius,
                 const double& max_interatomic_distance)
    : __score(nullptr), __ligand_idatm_types(nullptr) {
    try {
        __identify_gridpoints(centroids, grid, grid_spacing, dist_cutoff,
                              excluded_radius, max_interatomic_distance);
    } catch (...) {
        log_error << "FAILURE: scoreless constructor of Gpoints\n";
        throw;
    }
}

const Gpoints::Gpoint& Gpoints::get_center_point() const {
    geometry::Point center(0, 0, 0);
    double min_d = HUGE_VAL;
    const Gpoints::Gpoint* center_point = nullptr;
    for (auto& p : get_gridpoints0()) {
        const double d = p.crd().distance(center);
        if (d < min_d) {
            min_d = d;
            center_point = &p;
        }
    }
    if (!center_point) throw Error("die : central point cannot be found");
    return *center_point;
}

void Gpoints::__identify_gridpoints(const centro::Centroids& centroids,
                                    const molib::Atom::Grid& grid,
                                    const double& grid_spacing,
                                    const int& dist_cutoff,
                                    const double& excluded_radius,
                                    const double& max_interatomic_distance) {
#ifndef NDEGUG
    Benchmark bench;
#endif
    if (centroids.empty()) throw Error("die : there are no centroids");

    // find the absolute minimium and maximum coordinates of all centroids
    // find the absolute minimium and maximum coordinates of all centroids
    for (auto& kv : centroids) {
        geometry::Coordinate min(10000, 10000, 10000);
        geometry::Coordinate max(-10000, -10000, -10000);
        const int bsite_id = kv.first;
        for (auto& centroid : kv.second) {
            geometry::Coordinate min2 =
                centroid.get_centroid() - ceil(centroid.get_radial_check());
            geometry::Coordinate max2 =
                centroid.get_centroid() + ceil(centroid.get_radial_check());
            if (min2.x() < min.x()) min.set_x(min2.x());
            if (min2.y() < min.y()) min.set_y(min2.y());
            if (min2.z() < min.z()) min.set_z(min2.z());
            if (max2.x() > max.x()) max.set_x(max2.x());
            if (max2.y() > max.y()) max.set_y(max2.y());
            if (max2.z() > max.z()) max.set_z(max2.z());
        }
        dbgmsg("min point = " << min.pdb());
        dbgmsg("max point = " << max.pdb());
        const int total_gridpoints = 3 *
                                     ceil((max.x() - min.x()) / grid_spacing) *
                                     ceil((max.y() - min.y()) / grid_spacing) *
                                     ceil((max.z() - min.z()) / grid_spacing);
        log_note << "approximately " << total_gridpoints
                 << " gridpoints to evaluate\n\n";
        int model_number = 1;
        int points_kept = 0;
        int gridpoint_counter = 0;
        const double r = grid_spacing / 2;
        const double max_d = min.distance(max);  // distance between min and max
        const int last_column = ceil(max_d / r);
        const int last_row = ceil(max_d / (sqrt(3) * r));
        const int last_layer = ceil(max_d / (2 * r * sqrt(6) / 3));

        // initialize mapping between gridpoints and discretized 3D space
        __gmap[bsite_id].init(last_column + 1, last_row + 1, last_layer + 1);
        // dbgmsg("gmap szi = " << __gmap.szi << " szj = " << __gmap.szj << "
        // szk = " << __gmap.szk);
        geometry::Coordinate eval;
        for (int column = 0; column <= last_column; column++) {
            int even_column = (column % 2 == 0) ? 1 : 0;  // 1 if odd, 0 if even
            for (int row = 0; row <= last_row; row++) {
                int even_row = (row % 2 == 0) ? 1 : 0;  // 1 if odd, 0 if even
                for (int layer = 0; layer <= last_layer; layer++) {
                    int even_layer =
                        (layer % 2 == 0) ? 1 : 0;  // 1 if odd, 0 if even
                    if ((even_column == 0 && even_row == 0) ||
                        (even_column == 1 && even_row == 1)) {
                        if (even_layer == 1) {
                            eval.set_x(min.x() + column * r);
                            eval.set_y(min.y() + sqrt(3) * r * row);
                            eval.set_z(min.z() + layer * 2 * r * sqrt(6) / 3);
                        } else {
                            eval.set_x(min.x() + r + column * r);
                            eval.set_y(min.y() + r / sqrt(3) +
                                       sqrt(3) * r * row);
                            eval.set_z(min.z() + layer * 2 * r * sqrt(6) / 3);
                        }

                        // mark that this point is not accessible
                        //~ __gmap.data[column][row][layer] = false;
                        //~ assert(column < __gmap.szi);
                        //~ assert(row < __gmap.szj);
                        //~ assert(layer < __gmap.szk);
                        __gmap[bsite_id].data[column][row][layer] = nullptr;

                        int okay_min = 1;
                        int okay_max = 1;
                        double closest = 10000.0;
                        // if the point is within the radial_check of ANY of the
                        // centroids...
                        bool is_within_r_of_c = false;
                        for (auto& centroid : kv.second) {
                            if (eval.distance(centroid.get_centroid()) <=
                                centroid.get_radial_check()) {
                                is_within_r_of_c = true;
                                break;
                            }
                        }
                        if (is_within_r_of_c) {
                            for (molib::Atom* a :
                                 grid.get_neighbors(eval, dist_cutoff)) {
                                molib::Atom& atom = *a;
                                //~ dbgmsg("before getting atom radius");
                                const double vdW = atom.radius();
                                //~ dbgmsg("vdW = " << vdW);
                                const double eval_dist =
                                    excluded_radius + 0.9 * vdW;
                                const double distance =
                                    atom.crd().distance(eval);
                                if (distance <= eval_dist) {
                                    okay_min = 0;
                                    dbgmsg("distance = " << distance
                                                         << " eval_dist = "
                                                         << eval_dist);
                                    break;
                                } else {
                                    okay_min = 1;
                                    if (distance < closest) closest = distance;
                                }
                            }

                            if (closest > max_interatomic_distance)
                                okay_max = 0;
                            if (okay_min * okay_max > 0) {
                                dbgmsg("before adding to __gridpoints");
                                __gridpoints[bsite_id].push_back(Gpoint{
                                    eval, IJK{column, row, layer},
                                    __score
                                        ? __score->compute_energy(
                                              grid, eval, *__ligand_idatm_types)
                                        : Array1d<double>()});
                                dbgmsg("really out");
                                //~ __gmap.data[column][row][layer] = true;
                                model_number++;
                                points_kept++;
                            }
                        }
                    }
#ifndef NDEBUG
                    const int mod = gridpoint_counter % 10000;
                    if (mod == 0) {
                        double wall_secs = bench.seconds_from_start();
                        dbgmsg("Processing gridpoint "
                               << gridpoint_counter << " of approximately "
                               << total_gridpoints << " (took " << wall_secs
                               << " seconds)");
                        bench.reset();
                    }
#endif
                    gridpoint_counter++;
                }
                dbgmsg("column = " << column);
            }
        }
        // the last ones that did not get to the next mod==0
        log_note << points_kept << " points kept out of " << gridpoint_counter
                 << " total gridpoints\n";
    }
    // initialize gmap data here, because push_back can invalidate pointers...
    dbgmsg("initializing gmap");
    for (auto& kv : __gridpoints) {
        const int bsite_id = kv.first;
        dbgmsg("binding site = " << bsite_id
                                 << " number of points = " << kv.second.size());
        for (auto& gpoint : kv.second) {
            __gmap[bsite_id].data[gpoint.ijk().i][gpoint.ijk().j]
                                 [gpoint.ijk().k] = &gpoint;
        }
    }
}

void Gpoints::__identify_gridpoints(const double& grid_spacing,
                                    const double& radial_check) {
#ifndef NDEBUG
    Benchmark bench;
#endif
    // find the minimium and maximum coordinates
    geometry::Point center(0, 0, 0);

    geometry::Coordinate min = center - ceil(radial_check);
    geometry::Coordinate max = center + ceil(radial_check);
    dbgmsg("min point = " << min.pdb());
    dbgmsg("max point = " << max.pdb());
    const int total_gridpoints = 3 * ceil((max.x() - min.x()) / grid_spacing) *
                                 ceil((max.y() - min.y()) / grid_spacing) *
                                 ceil((max.z() - min.z()) / grid_spacing);
    log_note << "approximately " << total_gridpoints
             << " gridpoints to evaluate\n\n";
    int points_kept = 0;
    int gridpoint_counter = 0;
    const double r = grid_spacing / 2;
    const double max_d = min.distance(max);  // distance between min and max
    const int last_column = ceil(max_d / r);
    const int last_row = ceil(max_d / (sqrt(3) * r));
    const int last_layer = ceil(max_d / (2 * r * sqrt(6) / 3));

    geometry::Coordinate eval;
    for (int column = 0; column <= last_column; column++) {
        int even_column = (column % 2 == 0) ? 1 : 0;  // 1 if odd, 0 if even
        for (int row = 0; row <= last_row; row++) {
            int even_row = (row % 2 == 0) ? 1 : 0;  // 1 if odd, 0 if even
            for (int layer = 0; layer <= last_layer; layer++) {
                int even_layer =
                    (layer % 2 == 0) ? 1 : 0;  // 1 if odd, 0 if even
                if ((even_column == 0 && even_row == 0) ||
                    (even_column == 1 && even_row == 1)) {
                    if (even_layer == 1) {
                        eval.set_x(min.x() + column * r);
                        eval.set_y(min.y() + sqrt(3) * r * row);
                        eval.set_z(min.z() + layer * 2 * r * sqrt(6) / 3);
                    } else {
                        eval.set_x(min.x() + r + column * r);
                        eval.set_y(min.y() + r / sqrt(3) + sqrt(3) * r * row);
                        eval.set_z(min.z() + layer * 2 * r * sqrt(6) / 3);
                    }

                    // if the point is within the radial_check of ANY of the
                    // centroids...
                    if (eval.distance(center) <= radial_check) {
                        // Coordinare, IJK, Energy
                        __gridpoints[0].push_back(
                            Gpoint{eval, IJK{column, row, layer},
                                   Array1d<double>(0.0)});
                        dbgmsg("lastcolumn = "
                               << last_column << " lastrow = " << last_row
                               << " lastlayer = " << last_layer);
                        dbgmsg("column = " << column << " row = " << row
                                           << " layer = " << layer);
                        dbgmsg("gridpoint0 = " << __gridpoints[0].back().ijk());
                        points_kept++;
                    }
                }
#ifndef NDEBUG
                const int mod = gridpoint_counter % 10000;
                if (mod == 0) {
                    double wall_secs = bench.seconds_from_start();
                    dbgmsg("Processing gridpoint "
                           << gridpoint_counter << " of approximately "
                           << total_gridpoints << " (took " << wall_secs
                           << " seconds)");
                    bench.reset();
                }
#endif
                gridpoint_counter++;
            }
            dbgmsg("column = " << column);
        }
    }

    // center the ijk coordinates of grid points around the center point (0,0,0)
    Gpoint cp =
        this->get_center_point();  // here we copy by value intentionally
    dbgmsg("center point = " << cp.ijk());
    for (auto& point : this->get_gridpoints0()) {
        dbgmsg("point = " << point.ijk() << " address = " << &point);
        point.ijk() = point.ijk() - cp.ijk();
        dbgmsg("centered point = " << point.ijk());
    }

    // the last ones that did not get to the next mod==0
    log_note << points_kept << " points kept out of " << gridpoint_counter
             << " total gridpoints\n";
}
}
}
