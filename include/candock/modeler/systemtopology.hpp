#ifndef SYSTEMTOPOLOGY_H
#define SYSTEMTOPOLOGY_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "candock/geometry/geometry.hpp"
#include "candock/helper/debug.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/modeler/topology.hpp"
#include "openmm/Vec3.h"
#include "openmm/CustomNonbondedForce.h"

namespace OpenMM
{
class System;
class Integrator;
class Context;
class HarmonicAngleForce;
class HarmonicBondForce;
class PeriodicTorsionForce;
};

namespace candock {

namespace molib
{
class Atom;
class Molecule;
};

namespace OMMIface
{
struct ForceField;

class SystemTopology
{
      public:
        enum integrator_type
        {
                none,
                verlet,
                langevin,
                brownian,
        };

      private:
        OpenMM::System *system;
        OpenMM::Integrator *integrator;
        OpenMM::Context *context;
        OpenMM::CustomNonbondedForce *forcefield;

        integrator_type __integrator_used;
        int __thermostat_idx;

        OpenMM::HarmonicBondForce *bondStretch;
        OpenMM::HarmonicAngleForce *bondBend;
        OpenMM::PeriodicTorsionForce *bondTorsion;

        const ForceField *__ffield;

        int __kbforce_idx;
        std::vector<bool> masked;
        std::vector<double> masses;

        class AtomPoint
        {
              private:
                const geometry::Point __crd;
                molib::Atom &__atom;

              public:
                AtomPoint(const geometry::Point &crd, molib::Atom &atom) : __crd(crd), __atom(atom) {}
                const geometry::Point &crd() const
                {
                        return __crd;
                }
                molib::Atom &get_atom()
                {
                        return __atom;
                }
                void distance(double) const {} // just dummy : needed by grid

                typedef std::vector<std::unique_ptr<AtomPoint>> UPVec;
                typedef std::vector<AtomPoint *> PVec;
                typedef candock::molib::Grid<AtomPoint> Grid;
        };

        struct ForceData
        {
                int force_idx, idx1, idx2, idx3, idx4;
                double length, angle;
                int periodicity;
                double phase, k;
        };
        std::vector<std::vector<ForceData>> bondStretchData, bondBendData, bondTorsionData;

        void retype_amber_protein_atom_to_gaff(const molib::Atom &atom, int &type);

      public:
        SystemTopology() : system(nullptr), integrator(nullptr), context(nullptr),
                           __integrator_used(integrator_type::none), __thermostat_idx(-1), __kbforce_idx(-1) {}
        ~SystemTopology();
        static void loadPlugins(const std::string &extra_dir = "");
        void mask(Topology &topology, const molib::Atom::Vec &atoms);
        void unmask(Topology &topology, const molib::Atom::Vec &atoms);

        void mask_forces(const int atom_idx, const std::set<int> &substruct);
        void unmask_forces(const int atom_idx, const std::set<int> &substruct);

        void init_integrator(SystemTopology::integrator_type type,
                             const double step_size_in_ps,
                             const double temperature_in_kevin,
                             const double friction_in_per_ps);

        void init_platform(const std::string& platform, 
                           const std::string& precision,
                           const std::string& accelerators);

        void init_particles(Topology &topology);
        void init_physics_based_force(Topology &topology);
        void init_knowledge_based_force(Topology &topology);
        void init_bonded(Topology &topology, const bool use_constraints);
        void init_positions(const geometry::Point::Vec &crds);

        void update_thermostat(const double temperature_in_kevin,
                               const double collision_frequency);

        std::vector<OpenMM::Vec3> get_positions_in_nm();
        std::vector<OpenMM::Vec3> get_forces();
        double get_potential_energy();
        void minimize(const double tolerance, const int max_iterations);
        void dynamics(const int steps);
        void set_forcefield(const ForceField &ffield)
        {
                __ffield = &ffield;
        }
};
};

}

#endif
