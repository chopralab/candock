/* This is centroids.cpp and is part of CANDOCK
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

#include "candock/centro/centroids.hpp"
#include "candock/cluster/greedy.hpp"

#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/benchmark.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"

using namespace std;

namespace candock {
namespace centro {

using namespace statchem;

/* Centroid stuff
 *
 */
vector<Centroid> split_binding_site(
    const molib::Molecules& binding_site_ligands,
    const double centro_clus_rad) {
    vector<Centroid> centroids;
    geometry::Point::Vec crds = binding_site_ligands.get_crds();
    geometry::Point::Vec clustered =
        cluster::Cluster::greedy(crds, centro_clus_rad);
    for (auto& point : clustered) {
        // find closest point and calculate distance to this point
        geometry::Point closest;
        double min_dist = HUGE_VAL;
        for (auto& point2 : clustered) {
            const double d = point.distance(point2);
            if (d > 0 && d < min_dist) {
                min_dist = d;
                closest = point2;
            }
        }
        centroids.push_back(Centroid(
            point, point.distance(closest) + 1.0));  // add 1.0 A tolerance
        dbgmsg("ATOM      1  X   GEO X   1    "
               << point.pdb() << "  1.00  0.00           N  ");
    }
    return centroids;
}

Centroids set_centroids(
    const genlig::BindingSiteClusters& binding_site_clusters,
    const double centro_clus_rad) {
    Centroids centroids;
    if (binding_site_clusters.empty())
        throw Error(
            "die : no binding sites could be predicted for this protein - "
            "define its binding site(s) using the centroid option");
    for (auto& kv : binding_site_clusters) {
        const int bsite_id = kv.first;
        const molib::Molecules& binding_site_ligands = kv.second;
        centroids[bsite_id] =
            split_binding_site(binding_site_ligands, centro_clus_rad);
        dbgmsg("to better capture the shape of the binding site "
               << " it has been split into " << centroids[bsite_id].size()
               << " centroids");
    }
    dbgmsg("setting centroids from ProBiS predicted binding sites : "
           << endl
           << centroids);
    return centroids;
}
Centroids set_centroids(const string& centroid_file, const int num_bsites) {
    Centroids centroids;
    vector<string> data;
    fileio::read_file(centroid_file, data);
    for (string& line : data) {
        stringstream ss(line);
        int bsite_id;
        double x, y, z, rc;
        ss >> bsite_id >> x >> y >> z >> rc;
        if (bsite_id <= num_bsites) {
            centroids[bsite_id].push_back(
                Centroid(geometry::Coordinate(x, y, z), rc));
        }
    }
    if (centroids.empty())
        throw Error("die: could not find centroid in centroid file " +
                    centroid_file + "\n");
    dbgmsg("setting centroids from file : " << endl << centroids);
    return centroids;
}

ostream& operator<<(ostream& os, const centro::Centroids& centroids) {
    for (auto& kv : centroids) {
        const int bsite_id = kv.first;
        for (auto& centroid : kv.second) {
            os << bsite_id << " " << centroid.get_centroid().simple() << " "
               << fixed << setprecision(3) << centroid.get_radial_check()
               << endl;
        }
    }
    return os;
}
}
}
