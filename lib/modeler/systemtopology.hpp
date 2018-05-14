#ifndef SYSTEMTOPOLOGY_H
#define SYSTEMTOPOLOGY_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/geom3d.hpp"
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "molib/molecule.hpp"
#include "modeler/topology.hpp"
#include "openmm/Vec3.h"
#include <openmm/CustomNonbondedForce.h>

using namespace std;

namespace Molib
{
class Atom;
class Molecule;
};

namespace OpenMM
{
class System;
class Integrator;
class Context;
class HarmonicAngleForce;
class HarmonicBondForce;
class PeriodicTorsionForce;
};

namespace KBPlugin
{
class KBForce;
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

        //KBPlugin::KBForce* __kbforce;
        int __kbforce_idx;
        vector<bool> masked;
        vector<double> masses;

        class AtomPoint
        {
              private:
                const Geom3D::Point __crd;
                Molib::Atom &__atom;

              public:
                AtomPoint(const Geom3D::Point &crd, Molib::Atom &atom) : __crd(crd), __atom(atom) {}
                const Geom3D::Point &crd() const
                {
                        return __crd;
                }
                Molib::Atom &get_atom()
                {
                        return __atom;
                }
                void distance(double) const {} // just dummy : needed by grid

                typedef vector<unique_ptr<AtomPoint>> UPVec;
                typedef vector<AtomPoint *> PVec;
                typedef ::Grid<AtomPoint> Grid;
        };

        struct ForceData
        {
                int force_idx, idx1, idx2, idx3, idx4;
                double length, angle;
                int periodicity;
                double phase, k;
        };
        vector<vector<ForceData>> bondStretchData, bondBendData, bondTorsionData;

        void retype_amber_protein_atom_to_gaff(const Molib::Atom &atom, int &type);

      public:
        SystemTopology() : system(nullptr), integrator(nullptr), context(nullptr),
                           __integrator_used(integrator_type::none), __thermostat_idx(-1), __kbforce_idx(-1) {}
        ~SystemTopology();
        static void loadPlugins(const std::string &extra_dir = "");
        void mask(Topology &topology, const Molib::Atom::Vec &atoms);
        void unmask(Topology &topology, const Molib::Atom::Vec &atoms);

        void mask_forces(const int atom_idx, const set<int> &substruct);
        void unmask_forces(const int atom_idx, const set<int> &substruct);

        void init_integrator(SystemTopology::integrator_type type,
                             const double step_size_in_ps,
                             const double temperature_in_kevin,
                             const double friction_in_per_ps);

        void init_particles(Topology &topology);
        void init_physics_based_force(Topology &topology);
        void init_knowledge_based_force(Topology &topology);
        void init_bonded(Topology &topology, const bool use_constraints);
        void init_positions(const Geom3D::Point::Vec &crds);

        void update_thermostat(const double temperature_in_kevin,
                               const double collision_frequency);

        vector<OpenMM::Vec3> get_positions_in_nm();
        vector<OpenMM::Vec3> get_forces();
        double get_potential_energy();
        void minimize(const double tolerance, const double max_iterations);
        void dynamics(const int steps);
        void set_forcefield(const ForceField &ffield)
        {
                __ffield = &ffield;
        }
};
};
#endif
