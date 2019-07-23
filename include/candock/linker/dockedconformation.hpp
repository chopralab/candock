/* This is dockedconformation.hpp and is part of CANDOCK
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

#ifndef DOCKEDCONFORMATION_H
#define DOCKEDCONFORMATION_H
#include <memory>
#include "statchem/molib/molecule.hpp"

namespace candock {
namespace linker {

using namespace statchem;

class DockedConformation {
   public:
    typedef std::vector<DockedConformation> Vec;

   private:
    std::unique_ptr<molib::Molecule> __ligand;
    std::unique_ptr<molib::Molecule> __receptor;
    double __energy;
    double __potential_energy;
    size_t __max_clq_identity;

   public:
    DockedConformation()
        : __ligand(nullptr),
          __receptor(nullptr),
          __energy(0),
          __potential_energy(0),
          __max_clq_identity(0) {}

    DockedConformation(molib::Molecule ligand, molib::Molecule receptor,
                       double energy, double pot, size_t id)
        : __ligand(
              std::unique_ptr<molib::Molecule>(new molib::Molecule(ligand))),
          __receptor(
              std::unique_ptr<molib::Molecule>(new molib::Molecule(receptor))),
          __energy(energy),
          __potential_energy(pot),
          __max_clq_identity(id) {}

    molib::Molecule& get_ligand() { return *__ligand; }

    const molib::Molecule& get_ligand() const { return *__ligand; }

    molib::Molecule& get_receptor() { return *__receptor; }

    const molib::Molecule& get_receptor() const { return *__receptor; }

    double get_energy() const { return __energy; }

    double get_potential_energy() const { return __potential_energy; }

    size_t get_max_clq_identity() const { return __max_clq_identity; }

    void free_memory() {
        __ligand.reset(nullptr);
        __receptor.reset(nullptr);
    }

    friend bool operator<(const DockedConformation& other,
                          const DockedConformation& other2) {
        return other.get_energy() < other2.get_energy();
    }

    static void sort(DockedConformation::Vec& v);

    friend std::ostream& operator<<(std::ostream& os,
                                    const DockedConformation& conf);
};
};
}

#endif
