#include "modeler.hpp"
#include "forcefield.hpp"
#include "systemtopology.hpp"
#include "topology.hpp"
#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "helper/error.hpp"
#include "score/score.hpp"
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

#include "fileout/fileout.hpp"

using namespace std;

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


Modeler::Modeler(const ForceField &ffield, const string &fftype, double dist_cutoff,
          double tolerance, int max_iterations, int update_freq, double position_tolerance,
          bool use_constraints, double step_size_in_fs, double temperature, double friction)
      : __ffield(&ffield), __fftype(fftype), __dist_cutoff_in_nm(dist_cutoff * OpenMM::NmPerAngstrom),
        __tolerance(tolerance), __max_iterations(max_iterations), __update_freq(update_freq),
        __position_tolerance_in_nm(position_tolerance * OpenMM::NmPerAngstrom),
        __use_constraints(use_constraints), __step_size_in_ps(step_size_in_fs * OpenMM::PsPerFs),
        __temperature(temperature), __friction(friction), __run_dyanmics(false)
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

        void Modeler::add_crds(const Molib::Atom::Vec &atoms, const Geom3D::Point::Vec &crds) {
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
                        __positions[idx] = Geom3D::Point(x, y, z); 
                }
        }

/**
         * Changes coordinates of atoms
         */

Geom3D::Point::Vec Modeler::get_state(const Molib::Atom::Vec &atoms)
{

        const vector<OpenMM::Vec3> &positions_in_nm = __system_topology.get_positions_in_nm();
        Geom3D::Point::Vec crds;
        for (size_t i = 0; i < atoms.size(); ++i)
        {
                int idx = __topology.get_index(*atoms[i]);
                crds.push_back(Geom3D::Point(
                    positions_in_nm[idx][0] * OpenMM::AngstromsPerNm,
                    positions_in_nm[idx][1] * OpenMM::AngstromsPerNm,
                    positions_in_nm[idx][2] * OpenMM::AngstromsPerNm));
        }

        return crds;
}

void Modeler::minimize_state()
{
        __run_dyanmics = false;
        if (__fftype == "kb")
                knowledge_based_calculation();
        else if (__fftype == "phy")
                physical_calculation();
}

void Modeler::dynamics()
{
        __run_dyanmics = true;
        log_step << "Running " << __step_size_in_ps * __dynamics_steps << "ps simulation at " << __temperature << "K" << endl;
        if (__fftype == "kb")
                knowledge_based_calculation();
        else if (__fftype == "phy")
                physical_calculation();
}

void Modeler::physical_calculation()
{
        Benchmark bench;

        if (!__run_dyanmics)
        {
                log_step << "Doing energy minimization using physical forcefield" << endl;
                __system_topology.minimize(__tolerance, __max_iterations);
                log_benchmark << "time to minimize took " << bench.seconds_from_start()
                              << " wallclock seconds"
                              << "\n";
        }
        else
        {
                log_step << "Doing energy minimization using physical forcefield" << endl;
                __system_topology.dynamics(__dynamics_steps);
                log_benchmark << "time to minimize took " << bench.seconds_from_start()
                              << " wallclock seconds"
                              << "\n";
        }
}

void Modeler::knowledge_based_calculation()
{
        Benchmark bench;

        if (!__run_dyanmics)
        {
                __system_topology.minimize(__tolerance, __update_freq);
        }
        else
        {
                __system_topology.dynamics(__dynamics_steps);
        }

        log_benchmark << "Minimization took "
                      << bench.seconds_from_start() << " wallclock seconds"
                      << "\n";
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

void Modeler::init_openmm(SystemTopology::integrator_type type, string platform)
{
        __system_topology.set_forcefield(*__ffield);
        __system_topology.init_particles(__topology);
        __system_topology.init_bonded(__topology, __use_constraints);

        if (__fftype == "kb")
        {
                __system_topology.init_knowledge_based_force(__topology);
        }
        else if (__fftype == "phy")
        {
                __system_topology.init_physics_based_force(__topology);
        }
        else if (__fftype == "none")
        {
                //Do nothing
        }
        else
        {
                throw Error("die : unsupported forcefield");
        }

        __system_topology.init_integrator(type, __step_size_in_ps, __temperature, __friction, platform);
}

double Modeler::potential_energy()
{
        return __system_topology.get_potential_energy();
}
};
