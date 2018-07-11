#include "candock/modeler/topology.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <regex>
#include "candock/helper/debug.hpp"
#include "candock/helper/error.hpp"
#include "candock/helper/inout.hpp"
#include "candock/modeler/forcefield.hpp"
#include "candock/molib/molecule.hpp"
using namespace std;

namespace candock {
namespace OMMIface {
ostream& operator<<(ostream& stream,
                    const Topology::BondedExclusions& bonded_exclusions) {
    for (auto& atom_pair : bonded_exclusions) {
        auto& atom1 = *atom_pair.first;
        auto& atom2 = *atom_pair.second;
        stream << "bonded exclusions are atoms " << atom1.atom_number()
               << " and " << atom2.atom_number() << endl;
    }
    return stream;
}
ostream& operator<<(ostream& stream, const Topology::Bonds& bonds) {
    for (auto& atom_pair : bonds) {
        auto& atom1 = *atom_pair.first;
        auto& atom2 = *atom_pair.second;
        stream << "bond between atoms " << atom1.atom_number() << " and "
               << atom2.atom_number() << endl;
    }
    return stream;
}
ostream& operator<<(ostream& stream, const Topology::Angles& angles) {
    for (auto& angle : angles) {
        auto& atom1 = *get<0>(angle);
        auto& atom2 = *get<1>(angle);
        auto& atom3 = *get<2>(angle);
        stream << "angle between atoms " << atom1.atom_number() << " and "
               << atom2.atom_number() << " and " << atom3.atom_number() << endl;
    }
    return stream;
}
ostream& operator<<(ostream& stream, const Topology::Dihedrals& dihedrals) {
    for (auto& dihedral : dihedrals) {
        auto& atom1 = *get<0>(dihedral);
        auto& atom2 = *get<1>(dihedral);
        auto& atom3 = *get<2>(dihedral);
        auto& atom4 = *get<3>(dihedral);
        stream << "dihedral between atoms " << atom1.atom_number() << " and "
               << atom2.atom_number() << " and " << atom3.atom_number()
               << " and " << atom4.atom_number() << endl;
    }
    return stream;
}

int Topology::get_index(const molib::Atom& atom) const {
    if (!atom_to_index.count(&atom)) {
        log_error << "Problem with atom index: " << atom << endl;
    }

    return atom_to_index.at(&atom);
}

int Topology::get_type(const molib::Atom& atom) const {
    if (!atom_to_type.count(&atom)) {
        log_error << "Problem with atom type: " << atom << endl;
    }

    return atom_to_type.at(&atom);
}

Topology& Topology::add_topology(const molib::Atom::Vec& atoms,
                                 const ForceField& ffield) {
    int sz = this->atoms.size();
    this->atoms.insert(this->atoms.end(), atoms.begin(), atoms.end());

    // residue topology
    for (auto& patom : atoms) {
        molib::Atom& atom = *patom;
        const molib::Residue& residue = atom.br();

        if (!ffield.residue_topology.count(residue.resn()))
            throw Error("die : cannot find topology for residue " +
                        residue.resn());
        dbgmsg("residue topology for residue " << residue.resn());

        const ForceField::ResidueTopology& rtop =
            ffield.residue_topology.at(residue.resn());

        if (!rtop.atom.count(atom.atom_name()))
            throw Error("die : cannot find topology for atom " +
                        atom.atom_name() + " in residue " + residue.resn());

        dbgmsg("crd = " << atom.crd()
                        << " type = " << rtop.atom.at(atom.atom_name()));
        this->atom_to_type.insert({&atom, rtop.atom.at(atom.atom_name())});
        this->atom_to_index.insert({&atom, sz++});
    }

    BondedExclusions visited_bonds;
    set<tuple<molib::Atom*, molib::Atom*, molib::Atom*>> visited_angles;
    set<tuple<molib::Atom *, molib::Atom *, molib::Atom *, molib::Atom *>>
        visited_dihedrals, visited_impropers;

    // set the bonds, angles, dihedrals
    for (auto& patom1 : atoms) {
        auto& atom1 = *patom1;
        for (auto& atom2 : atom1) {
            if (!visited_bonds.count({&atom1, &atom2})) {
                visited_bonds.insert({&atom2, &atom1});
                this->bonds.push_back({&atom1, &atom2});
                dbgmsg("info.bond = " << atom1.atom_number() << " "
                                      << atom2.atom_number());
            }
            for (auto& atom3 : atom2) {
                if (&atom3 != &atom1) {
                    if (!visited_angles.count(
                            make_tuple(&atom1, &atom2, &atom3))) {
                        visited_angles.insert(
                            make_tuple(&atom3, &atom2, &atom1));
                        this->angles.push_back(
                            make_tuple(&atom1, &atom2, &atom3));
                        dbgmsg("info.angle = " << atom1.atom_number() << " "
                                               << atom2.atom_number() << " "
                                               << atom3.atom_number());
                    }
                    for (auto& atom4 : atom3) {  // propers
                        if (&atom4 != &atom2 && &atom4 != &atom1) {
                            if (!visited_dihedrals.count(make_tuple(
                                    &atom1, &atom2, &atom3, &atom4))) {
                                visited_dihedrals.insert(
                                    make_tuple(&atom4, &atom3, &atom2, &atom1));
                                this->dihedrals.push_back(
                                    make_tuple(&atom1, &atom2, &atom3, &atom4));
                                dbgmsg("info.dihedral (proper) = "
                                       << atom1.atom_number() << " "
                                       << atom2.atom_number() << " "
                                       << atom3.atom_number() << " "
                                       << atom4.atom_number());
                            }
                        }
                    }
                    for (auto& atom4 : atom2) {  // impropers
                        if (&atom4 != &atom3 && &atom4 != &atom1) {
                            if (!(visited_impropers.count(make_tuple(
                                      &atom3, &atom1, &atom2, &atom4)) ||
                                  visited_impropers.count(make_tuple(
                                      &atom1, &atom3, &atom2, &atom4)) ||
                                  visited_impropers.count(make_tuple(
                                      &atom1, &atom4, &atom2, &atom3)) ||
                                  visited_impropers.count(make_tuple(
                                      &atom4, &atom1, &atom2, &atom3)) ||
                                  visited_impropers.count(make_tuple(
                                      &atom4, &atom3, &atom2, &atom1)) ||
                                  visited_impropers.count(make_tuple(
                                      &atom3, &atom4, &atom2, &atom1)))) {
                                visited_impropers.insert(
                                    make_tuple(&atom4, &atom3, &atom2, &atom1));
                                this->impropers.push_back(
                                    make_tuple(&atom3, &atom1, &atom2, &atom4));
                                dbgmsg(
                                    "info.improper (improper) (third atom is "
                                    "central) = "
                                    << atom3.atom_number() << " "
                                    << atom1.atom_number() << " "
                                    << atom2.atom_number() << " "
                                    << atom4.atom_number());
                            }
                        }
                    }
                }
            }
        }
    }

    // set bonded exclusions
    for (auto& patom1 : atoms) {
        molib::Atom& atom1 = *patom1;
        for (auto& atom2 : atom1) {
            this->bonded_exclusions.insert({&atom1, &atom2});
            for (auto& atom3 : atom2) {
                if (&atom3 != &atom1) {
                    this->bonded_exclusions.insert({&atom1, &atom3});
                    for (auto& atom4 : atom3) {
                        if (&atom4 != &atom1 && &atom4 != &atom2) {
                            this->bonded_exclusions.insert({&atom1, &atom4});
                        }
                    }
                }
            }
        }
    }

    dbgmsg("BONDS ARE : " << endl << this->bonds);
    dbgmsg("ANGLES ARE : " << endl << this->angles);
    dbgmsg("DIHEDRALS ARE : " << endl << this->dihedrals);
    dbgmsg("IMPROPERS ARE : " << endl << this->impropers);
    dbgmsg("BONDED EXCLUSIONS ARE : " << endl << this->bonded_exclusions);

    return *this;
}
};
}
