/* This is conformations.hpp and is part of CANDOCK
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

#ifndef CONFORMATIONS_H
#define CONFORMATIONS_H

#include "statchem/helper/array2d.hpp"
#include "candock/docker/gpoints.hpp"

namespace candock {

namespace docker {

class Conformations {
   public:
    typedef std::map<int, std::map<int, std::map<int, std::vector<int>>>>
        ConfMap;

   private:
    const molib::Molecule& __seed;

    std::vector<Gpoints::PGpointVec> __conf_vec;
    ConfMap __conf_map;
    Array2d<double> __rmsd_sq;

   public:
    Conformations(const molib::Molecule& seed, const Gpoints& gpoints,
                  const double& conf_spin, const int num_univec);

    std::vector<Gpoints::PGpointVec>& get_conformations() { return __conf_vec; }
    const std::vector<Gpoints::PGpointVec>& get_conformations() const {
        return __conf_vec;
    }
    ConfMap& get_confmap() { return __conf_map; }
    std::vector<int>& get_confs_at(const Gpoints::IJK& ijk) {
        return __conf_map[ijk.i][ijk.j][ijk.k];
    }
    double get_rmsd_sq(const int i, const int j) {
        return __rmsd_sq.data[i][j];
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const Conformations& conformations);
};
}
}

#endif
