/* This is dock.hpp and is part of CANDOCK
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

#ifndef DOCK_H
#define DOCK_H

#include "candock/docker/conformations.hpp"
#include "candock/docker/gpoints.hpp"
#include "statchem/molib/molecules.hpp"

namespace candock {
namespace docker {

using namespace statchem;

class Dock {
    class DockedConf {
       public:
        typedef std::vector<DockedConf> Vec;

       private:
        const Gpoints::Gpoint& __cavpoint;
        Gpoints::PGpointVec& __conf0;
        double __energy;
        size_t __i;
        size_t __bsite_id;

       public:
        DockedConf(const Gpoints::Gpoint& cavpoint, Gpoints::PGpointVec& conf0,
                   double energy, size_t i, size_t bsite_id)
            : __cavpoint(cavpoint),
              __conf0(conf0),
              __energy(energy),
              __i(i),
              __bsite_id(bsite_id) {}
        // geometry::Point& crd() { return __cavpoint.crd(); }
        const geometry::Point& crd() const { return __cavpoint.crd(); }
        void distance(const double) const {}  // dummy
        const Gpoints::Gpoint& get_cavpoint() const { return __cavpoint; }
        const Gpoints::PGpointVec& get_conf0() const { return __conf0; }
        double get_energy() const { return __energy; }
        size_t get_i() const { return __i; }
        size_t get_bsite_id() const { return __bsite_id; }
        struct by_energy {
            bool operator()(const DockedConf* lhs,
                            const DockedConf* rhs) const {
                return lhs->__energy < rhs->__energy;
            }
        };
        double compute_rmsd_sq(const DockedConf& other) const;
    };

    const Gpoints& __gpoints;
    Conformations& __conformations;
    const molib::Molecule& __seed;
    double __rmsd_tol_sq;
    double __rmsd_tol;

    molib::Molecules __docked;

#ifndef NDEBUG
    const score::Score& __score;
    const molib::Atom::Grid& __gridrec;
#endif
    DockedConf::Vec __dock();
    DockedConf::Vec __cluster(const DockedConf::Vec& confs);
    void __cluster_fast(const DockedConf::Vec& conformations,
                        DockedConf::Vec& reps);
    void __set_docked(const DockedConf::Vec& confs);

   public:
#ifndef NDEBUG
    Dock(const Gpoints& gpoints, Conformations& conformations,
         const molib::Molecule& seed, const score::Score& score,
         const molib::Atom::Grid& gridrec, const double rmsd_tol = 2.0)
        : __gpoints(gpoints),
          __conformations(conformations),
          __seed(seed),
          __rmsd_tol_sq(rmsd_tol * rmsd_tol),
          __rmsd_tol(rmsd_tol),
          __score(score),
          __gridrec(gridrec) {}
#else
    Dock(const Gpoints& gpoints, Conformations& conformations,
         molib::Molecule& seed, const double rmsd_tol = 2.0)
        : __gpoints(gpoints),
          __conformations(conformations),
          __seed(seed),
          __rmsd_tol_sq(rmsd_tol * rmsd_tol),
          __rmsd_tol(rmsd_tol) {}
#endif
    void run();
    molib::Molecules& get_docked() { return __docked; };
};
}
}

#endif
