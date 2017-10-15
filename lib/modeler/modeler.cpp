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

#include "fileout/fileout.hpp"

using namespace std;

namespace OMMIface {
        ostream& operator<<(ostream& os, const vector<OpenMM::Vec3>& positions)	{
                os << "SIZE OF POSITIONS = " << positions.size() << endl;
                os << "ELEMENTS :" << endl;
                for (auto &v : positions) 
                        os << setprecision(8) << fixed << v[0] << " " << v[1] << " " << v[2] << endl;
                return os;
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

        Geom3D::Point::Vec Modeler::get_state(const Molib::Atom::Vec &atoms) {

                const vector<OpenMM::Vec3>& positions_in_nm = __system_topology.get_positions_in_nm();
                Geom3D::Point::Vec crds;
                for (size_t i = 0; i < atoms.size(); ++i) {
                        int idx = __topology.get_index(*atoms[i]);
                        crds.push_back(Geom3D::Point(
                                positions_in_nm[idx][0] * OpenMM::AngstromsPerNm,
                                positions_in_nm[idx][1] * OpenMM::AngstromsPerNm,
                                positions_in_nm[idx][2] * OpenMM::AngstromsPerNm
                        ));
                }

                return crds;
        }

        void Modeler::minimize_state() {
                __run_dyanmics = false;
                if (__fftype == "kb")
                        knowledge_based_calculation();
                else if (__fftype == "phy")
                        physical_calculation();
        }

        void Modeler::dynamics() {
                __run_dyanmics = true;
                log_step << "Running " << __step_size_in_ps * __dynamics_steps << "ps simulation at " << __temperature << "K" << endl;
                if (__fftype == "kb")
                        knowledge_based_calculation();
                else if (__fftype == "phy")
                        physical_calculation();
        }

        void Modeler::physical_calculation() {
                Benchmark bench;

                if (! __run_dyanmics) {
                        log_step << "Doing energy minimization using physical forcefield" << endl;
                        __system_topology.minimize(__tolerance, __max_iterations);
                        log_benchmark << "time to minimize took " << bench.seconds_from_start() 
                                      << " wallclock seconds" << "\n";
                } else {
                        log_step << "Doing energy minimization using physical forcefield" << endl;
                        __system_topology.dynamics(__dynamics_steps);
                        log_benchmark << "time to minimize took " << bench.seconds_from_start() 
                                      << " wallclock seconds" << "\n";
                }

        }

        void Modeler::knowledge_based_calculation() {
                Benchmark bench;
                dbgmsg( "Doing energy minimization using knowledge-based forcefield");

                // for knowledge-based forcefield we implement a custom update nonbond
                // function

                int iter = 0;

                // do a brief relaxation of bonded forces initially (without non-bonded forces)
                __system_topology.clear_knowledge_based_force();
                __system_topology.minimize(__tolerance, 20);

                vector<OpenMM::Vec3> initial_positions = __system_topology.get_positions_in_nm();
                dbgmsg("initial_positions (after brief minimization) = " << initial_positions);

                // check if minimization failed
                if (std::isnan(initial_positions[0][0]))
                        throw MinimizationError("die : minimization failed (initial bonded relaxation)");

                while (iter < __max_iterations) {

                        dbgmsg("starting minimization step = " << iter);

                        dbgmsg("initial_positions = " << initial_positions);

                        __system_topology.update_knowledge_based_force(__topology, initial_positions, __dist_cutoff_in_nm);

                        if (! __run_dyanmics) {
                                __system_topology.minimize(__tolerance, __update_freq);
                        } else {
                                __system_topology.dynamics(__dynamics_steps);
                        }

                        const vector<OpenMM::Vec3>& minimized_positions = __system_topology.get_positions_in_nm();

                        // check if minimization failed
                        if (std::isnan(minimized_positions[0][0]))
                                throw MinimizationError("die : minimization failed (in loop)");

                        dbgmsg("minimized_positions = " << minimized_positions);

                        // check if positions have converged
                        double max_error = 0;
                        for (size_t i = 0; i < initial_positions.size(); ++i) {
                                OpenMM::Vec3 dif_pos = initial_positions[i] - minimized_positions[i];
                                const double error = dif_pos.dot(dif_pos);
                                    if (error > max_error)
                                        max_error = error; 
                        }

                        // stop if convergence reached
                        if (sqrt(max_error) < __position_tolerance_in_nm) {
                                dbgmsg("Convergence reached (position_tolerance) - minimization finished");
                                break;
                        }

                        iter += __update_freq;

                        initial_positions = minimized_positions;
                        dbgmsg("ending minimization step iter = " << iter);

                }

                __system_topology.clear_knowledge_based_force();
                log_benchmark << "Minimized in " << iter << " iterations, which took " 
                        << bench.seconds_from_start() << " wallclock seconds" << "\n";
        }

#ifndef NDEBUG
        void Modeler::minimize_knowledge_based(Molib::Molecule& ligand, Molib::Molecule& receptor, Score::Score& score) {
                Benchmark bench;
                dbgmsg( "Doing energy minimization of ligand " << ligand.name() << " using knowledge-based forcefield");

                // for knowledge-based forcefield we implement a custom update nonbond
                // function

                int iter = 0;

                // do a brief relaxation of bonded forces initially (without non-bonded forces)
                __system_topology.clear_knowledge_based_force();
                __system_topology.minimize(__tolerance, 20);

                vector<OpenMM::Vec3> initial_positions = __system_topology.get_positions_in_nm();
                dbgmsg("initial_positions (after brief minimization) = " << initial_positions);

                // check if minimization failed
                if (std::isnan(initial_positions[0][0]))
                        throw MinimizationError("die : minimization failed (initial bonded relaxation)");

                __system_topology.update_knowledge_based_force(__topology, initial_positions, __dist_cutoff_in_nm);

                while (iter < __max_iterations) {

                        dbgmsg("starting minimization step = " << iter);

                        dbgmsg("initial_positions = " << initial_positions);

                        // output frames during minimization
                        Molib::Molecule minimized_receptor(receptor, get_state(receptor.get_atoms()));
                        Molib::Molecule minimized_ligand(ligand, get_state(ligand.get_atoms()));

                        minimized_receptor.undo_mm_specific();

                        Molib::Atom::Grid gridrec(minimized_receptor.get_atoms());
                        const double energy = score.non_bonded_energy(gridrec, minimized_ligand);


                        std::stringstream ss;
                        Fileout::print_complex_pdb(ss, minimized_ligand, minimized_receptor, energy);
                        Inout::output_file(ss.str(), ligand.name() + "_frame_" + std::to_string(iter) + ".pdb");

                        __system_topology.minimize(__tolerance, __update_freq);

                        const vector<OpenMM::Vec3>& minimized_positions = __system_topology.get_positions_in_nm();

                        // check if minimization failed
                        if (std::isnan(minimized_positions[0][0]))
                                throw MinimizationError("die : minimization failed (in loop)");

                        dbgmsg("minimized_positions = " << minimized_positions);

                        const vector<OpenMM::Vec3>& forces = __system_topology.get_forces();
                        dbgmsg("forces after minimization = " << forces);

                        // check if positions have converged
                        double max_error = 0;
                        for (size_t i = 0; i < initial_positions.size(); ++i) {
                                OpenMM::Vec3 dif_pos = initial_positions[i] - minimized_positions[i];
                                const double error = dif_pos.dot(dif_pos);
                                if (error > max_error)
                                        max_error = error;
                        }

                        // stop if convergence reached
                        if (sqrt(max_error) < __position_tolerance_in_nm) {
                                dbgmsg("Convergence reached (position_tolerance) - minimization finished");
                                break;
                        }

                        // update knowledge-based nonbond list
                        __system_topology.update_knowledge_based_force(__topology, minimized_positions, __dist_cutoff_in_nm);

                        iter += __update_freq;

                        initial_positions = minimized_positions;
                        dbgmsg("ending minimization step iter = " << iter);

                }
                __system_topology.clear_knowledge_based_force();
                log_benchmark << "Minimized " << ligand.name() << " and " << receptor.name()
                        << " in " << iter << " iterations, which took " 
                        << bench.seconds_from_start() << " wallclock seconds" << "\n";
        }
#endif

        void Modeler::init_openmm_positions() {
#ifndef NDEBUG
                dbgmsg("init_openmm_positions");
                for (auto &point : __positions) dbgmsg(point);
#endif
                __system_topology.init_positions(__positions);
        }

        void Modeler::init_openmm(SystemTopology::integrator_type type) {
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
        }

        double Modeler::potential_energy() {
                return __system_topology.get_potential_energy();
        }

};
