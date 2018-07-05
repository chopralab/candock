#include "candock/modeler/modeler.hpp"
#include "candock/modeler/forcefield.hpp"
#include "candock/modeler/systemtopology.hpp"
#include "candock/modeler/topology.hpp"
#include "candock/helper/inout.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/error.hpp"
#include "candock/score/score.hpp"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <stdlib.h>
#include <time.h>

#include <openmm/Units.h>
#include <openmm/Vec3.h>

#include "candock/fileout/fileout.hpp"

using namespace std;

namespace candock {
namespace OMMIface
{
ostream &operator<<(ostream &os, const vector<OpenMM::Vec3> &positions)
{
        os << "SIZE OF POSITIONS = " << positions.size() << endl;
        os << "ELEMENTS :" << endl;
        for (auto &v : positions)
                os << setprecision(8) << fixed << v[0] << " " << v[1] << " " << v[2] << endl;
        return os;
}


Modeler::Modeler(
        const ForceField &ffield,
        const string &fftype,
        double tolerance,
        int max_iterations,
        bool use_constraints,
        double step_size_in_fs,
        double temperature,
        double friction)
      : __ffield(&ffield),
        __fftype(fftype),
        __tolerance(tolerance),
        __max_iterations(max_iterations),
        __use_constraints(use_constraints),
        __step_size_in_ps(step_size_in_fs * OpenMM::PsPerFs),
        __temperature(temperature),
        __friction(friction)
{
}

        void Modeler::mask(const Molib::Atom::Vec &atoms) {
                dbgmsg("Masking atoms " << atoms);
                __system_topology.mask(__topology, atoms);
        }

        void Modeler::unmask(const Molib::Atom::Vec &atoms) {
                dbgmsg("Unmasking atoms " << atoms);
                __system_topology.unmask(__topology, atoms);
        }

        void Modeler::add_topology(const Molib::Atom::Vec &atoms) {
                __topology.add_topology(atoms, *__ffield);
                __positions.resize(__topology.atoms.size()); // as many positions as there are atoms
        }

        void Modeler::add_crds(const Molib::Atom::Vec &atoms, const geometry::Point::Vec &crds) {
                for (size_t i = 0; i < atoms.size(); ++i) {
                        int idx = __topology.get_index(*atoms[i]);
                        __positions[idx] = crds[i];
                }
        }

        void Modeler::add_random_crds(const Molib::Atom::Vec &atoms) {
                srand(time(NULL));
                for (size_t i = 0; i < atoms.size(); ++i) {
                        int idx = __topology.get_index(*atoms[i]);
                        const double x = rand() % 100 + 1;
                        const double y = rand() % 100 + 1;
                        const double z = rand() % 100 + 1;
                        __positions[idx] = geometry::Point(x, y, z); 
                }
        }

/**
         * Changes coordinates of atoms
         */

geometry::Point::Vec Modeler::get_state(const Molib::Atom::Vec &atoms)
{

        const vector<OpenMM::Vec3> &positions_in_nm = __system_topology.get_positions_in_nm();
        geometry::Point::Vec crds;
        for (size_t i = 0; i < atoms.size(); ++i)
        {
                int idx = __topology.get_index(*atoms[i]);
                crds.push_back(geometry::Point(
                    positions_in_nm[idx][0] * OpenMM::AngstromsPerNm,
                    positions_in_nm[idx][1] * OpenMM::AngstromsPerNm,
                    positions_in_nm[idx][2] * OpenMM::AngstromsPerNm));
        }

        return crds;
}

void Modeler::minimize_state()
{
        Benchmark bench;
        __system_topology.minimize(__tolerance, __max_iterations);
        log_benchmark << "Minimization took "
                      << bench.seconds_from_start() << " wallclock seconds"
                      << "\n";
}

void Modeler::dynamics()
{
        log_step << "Running " << __step_size_in_ps * __dynamics_steps << "ps simulation at " << __temperature << "K" << endl;
        __system_topology.dynamics(__dynamics_steps);
}

void Modeler::init_openmm_positions()
{
#ifndef NDEBUG
        dbgmsg("init_openmm_positions");
        for (auto &point : __positions)
                dbgmsg(point);
#endif
        __system_topology.init_positions(__positions);
}

void Modeler::init_openmm(const std::string& platform, const std::string& precision, const std::string& accelerators, SystemTopology::integrator_type type)
{
        __system_topology.set_forcefield(*__ffield);
        __system_topology.init_particles(__topology);
        __system_topology.init_bonded(__topology, __use_constraints);

        if (__fftype == "kb") {
                __system_topology.init_knowledge_based_force(__topology);
        } else if (__fftype == "phy") {
                __system_topology.init_physics_based_force(__topology);
        } else if (__fftype == "none") {
                //Do nothing
        } else {
                throw Error("die : unsupported forcefield");
        }

        __system_topology.init_integrator(type, __step_size_in_ps, __temperature, __friction);
        __system_topology.init_platform(platform, precision, accelerators);
}

double Modeler::potential_energy()
{
        return __system_topology.get_potential_energy();
}
};
}
