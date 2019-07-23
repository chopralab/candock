/* This is greedy.hpp and is part of CANDOCK
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

#ifndef GREEDY_CLUSTER_H
#define GREEDY_CLUSTER_H

#include "statchem/geometry/geometry.hpp"
#include "statchem/score/score.hpp"
#include "candock/linker/partial.hpp"

namespace candock {
namespace cluster {

using namespace statchem;

class Cluster {
    template <typename T>
    class LinkedConf {
       public:
        typedef std::vector<std::unique_ptr<LinkedConf>> UPVec;

       private:
        T& __molecule;
        geometry::Point __crd;
        double __energy;

       public:
        LinkedConf(T& molecule, geometry::Point crd, double energy)
            : __molecule(molecule), __crd(crd), __energy(energy) {}
        geometry::Point& crd() { return __crd; }
        const geometry::Point& crd() const { return __crd; }
        void distance(const double) const {}  // dummy
        T& get_molecule() const { return __molecule; }
        double get_energy() const { return __energy; }
        struct by_energy {
            bool operator()(const LinkedConf* lhs,
                            const LinkedConf* rhs) const {
                return lhs->__energy < rhs->__energy;
            }
        };
    };
    friend std::ostream& operator<<(
        std::ostream& os,
        const std::set<const LinkedConf<molib::Molecule>*,
                       LinkedConf<molib::Molecule>::by_energy>& confs);
    friend std::ostream& operator<<(
        std::ostream& os,
        const std::set<const LinkedConf<linker::Partial>*,
                       LinkedConf<linker::Partial>::by_energy>& confs);

   public:
    static geometry::Point::Vec greedy(const geometry::Point::Vec& initial,
                                       const double clus_rad);
    static molib::Molecules greedy(const molib::Molecules& initial,
                                   const score::Score& score,
                                   molib::Atom::Grid& gridrec,
                                   const double clus_rad);
    static linker::Partial::Vec greedy(const linker::Partial::Vec& initial,
                                       const molib::Atom::Grid& gridrec,
                                       const double clus_rad);
};
};
}

#endif
