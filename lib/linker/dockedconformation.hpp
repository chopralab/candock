#ifndef DOCKEDCONFORMATION_H
#define DOCKEDCONFORMATION_H
#include "molib/molecule.hpp"
#include <memory>

namespace Linker {

        class DockedConformation {
        public:
                typedef vector<DockedConformation> Vec;
        private:
                unique_ptr<Molib::Molecule> __ligand;
                unique_ptr<Molib::Molecule> __receptor;
                double __energy;
                double __potential_energy;
                size_t __max_clq_identity;
        public:
                DockedConformation() : __ligand (nullptr), __receptor (nullptr), __energy (0), __potential_energy (0), __max_clq_identity(0) {}

                DockedConformation (Molib::Molecule ligand, Molib::Molecule receptor, double energy, double pot, size_t id)
                        : __ligand (unique_ptr<Molib::Molecule> (new Molib::Molecule (ligand))),
                          __receptor (unique_ptr<Molib::Molecule> (new Molib::Molecule (receptor))),
                          __energy (energy), __potential_energy (pot), __max_clq_identity(id) {}

                Molib::Molecule &get_ligand() {
                        return *__ligand;
                }

                const Molib::Molecule &get_ligand() const {
                        return *__ligand;
                }

                Molib::Molecule &get_receptor() {
                        return *__receptor;
                }

                const Molib::Molecule &get_receptor() const {
                        return *__receptor;
                }

                double get_energy() const {
                        return __energy;
                }

                double get_potential_energy() const {
                        return __potential_energy;
                }

                size_t get_max_clq_identity() const {
                        return __max_clq_identity;
                }

                friend bool operator< (const DockedConformation &other, const DockedConformation &other2) {
                        return other.get_energy() < other2.get_energy();
                }

                static void sort (DockedConformation::Vec &v);

                friend ostream &operator<< (ostream &os, const DockedConformation &conf);
        };

};

#endif
