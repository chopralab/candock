#include "modeler.hpp"
#include "forcefield.hpp"
#include "systemtopology.hpp"
#include "topology.hpp"
#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "helper/error.hpp"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
using namespace std;

namespace OMMIface {
	void Modeler::mask(const Molib::Atom::Vec &atoms) {
		__system_topology.mask(__topology, atoms);
	}

	void Modeler::unmask(const Molib::Atom::Vec &atoms) {
		__system_topology.unmask(__topology, atoms);
	}


	void Modeler::add_topology(Molib::Atom::Vec &atoms) {
		__positions.(__positions.size() + atoms.size()); // increase size to accomodate more atoms
		__topology.add_topology(atoms, *__ffield);
	}

	void Modeler::add_crds(Molib::Atom::Vec &atoms, Geom3D::Point::Vec &crds) {
		for (int i = 0; i < atoms.size(); ++i) {
			int idx = __topology.get_index(atoms[i]);
			__positions[idx] = crds[i];
		}
	}

	/**
	 * Changes coordinates of atoms
	 */
	Molib::Atom::Vec& Modeler::get_state(Molib::Atom::Vec &atoms) {

		int infoMask = OpenMM::State::Positions;
		const vector<OpenMM::Vec3>& positions_in_nm = __system_topology.get_positions_in_nm();

		for (int i = 0; i < atoms.size(); ++i) {
			int idx = __topology.get_index(atoms[i]);
			atoms[i].set_crd(Geom3D::Coordinate(
				positions_in_nm[idx][0] * OpenMM::AngstromsPerNm,
				positions_in_nm[idx][1] * OpenMM::AngstromsPerNm,
				positions_in_nm[idx][2] * OpenMM::AngstromsPerNm
				));
		}
		
		return atoms;
	}

	void Modeler::minimize_state() {
		if (__fftype == "kb")
			minimize_knowledge_based();
		else if (__fftype == "phy")
			minimize_physical();
	}
	
	void Modeler::minimize_physical() {
		Benchmark::reset();
		cout << "Doing energy minimization using physical forcefield" << endl;

		__system_topology.minimize(__tolerance, __max_iterations);

		cout << "time to minimize took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
	}
	
	void Modeler::minimize_knowledge_based() {
		Benchmark::reset();
		cout << "Doing energy minimization using knowledge-based forcefield" << endl;

		// for knowledge-based forcefield we implement a custom update nonbond
		// function
		
		int iter = 0;

		const OpenMM::vector<Vec3>& initial_positions = __system_topology.get_positions_in_nm();

		__system_topology.update_knowledge_based_force(__topology, initial_positions, __dist_cutoff_in_nm);

		while (iter < __max_iterations) {
			
			dbgmsg("starting minimization step = " << iter);

			__system_topology.minimize(__tolerance, __update_freq);

			const OpenMM::vector<Vec3>& minimized_positions = __system_topology.get_positions_in_nm();

			// check if positions have converged
			double max_error = 0;
			for (int i = 0; i < initial_positions.size(); ++i) {
				OpenMM::Vec3 dif_pos = initial_positions[i] - minimized_positions[i];
				 const double error = dif_pos.dot(dif_pos);
				 if (error > max_error)
					max_error = error;				 
			}

			// stop if convergence reached
			if (sqrt(max_error) < __position_tolerance_in_nm)
				break;

			// update knowledge-based nonbond list
			__system_topology.update_knowledge_based_force(__topology, minimized_positions, __dist_cutoff_in_nm);
				
			iter += update_freq;						

			initial_positions = minimized_positions;
			dbgmsg("ending minimization step iter = " << iter);

		}
		cout << "time to minimize took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
	}

	void Modeler::init_openmm() {
		__system_topology = SystemTopology::create_object(__ffield, 
			__positions, __topology, __step_size_in_fs, __use_constraints);
			
		__system_topology.init_integrator(__step_size_in_fs);
		__system_topology.init_positions(__positions);
		__system_topology.init_particles(__topology);
		__system_topology.init_bonded(__topology, __use_constraints);
		
		if (__fftype == "kb") {

			// do nothing since force will be initialized when minimization starts
			
		} else if (__fftype == "phy") {
			
			__system_topology.init_physics_based_force(__topology);
			
		} else {
			throw Error("die : unsupported forcefield");
		}
		
	}

};
