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
		for (int i = 0; i < atoms.size(); ++i) {
			int idx = __topology.get_index(*atoms[i]);
			__positions[idx] = crds[i];
		}
	}

	/**
	 * Changes coordinates of atoms
	 */

	Geom3D::Point::Vec Modeler::get_state(const Molib::Atom::Vec &atoms) {

		const vector<OpenMM::Vec3>& positions_in_nm = __system_topology.get_positions_in_nm();
		Geom3D::Point::Vec crds;
		for (int i = 0; i < atoms.size(); ++i) {
			int idx = __topology.get_index(*atoms[i]);
			crds.push_back(Geom3D::Point(
				positions_in_nm[idx][0] * OpenMM::AngstromsPerNm,
				positions_in_nm[idx][1] * OpenMM::AngstromsPerNm,
				positions_in_nm[idx][2] * OpenMM::AngstromsPerNm
			));
		}
		
		return crds;
	}


	void Modeler::minimize_state(const Molib::Molecule &ligand, const Molib::Molecule &receptor, const Molib::Score &score) {
		if (__fftype == "kb")
			minimize_knowledge_based(ligand, receptor, score);
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

	void Modeler::minimize_knowledge_based(const Molib::Molecule &ligand, const Molib::Molecule &receptor, const Molib::Score &score) {
		Benchmark::reset();
		cout << "Doing energy minimization of ligand " << ligand.name() << " using knowledge-based forcefield" << endl;

		// for knowledge-based forcefield we implement a custom update nonbond
		// function
		
		int iter = 0;

		// do a brief relaxation of bonded forces initially (without non-bonded forces)
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



#ifndef NDEBUG
			// output frames during minimization
			Molib::Molecule minimized_receptor(receptor, get_state(receptor.get_atoms()));
			Molib::Molecule minimized_ligand(ligand, get_state(ligand.get_atoms()));

			minimized_receptor.undo_mm_specific();

			Molib::Atom::Grid gridrec(minimized_receptor.get_atoms());
			const double energy = score.non_bonded_energy(gridrec, minimized_ligand);


			inout::output_file(Molib::Molecule::print_complex(minimized_ligand, minimized_receptor, energy), 
				ligand.name() + "_frame_" + help::to_string(iter) + ".pdb");
#endif

			__system_topology.minimize(__tolerance, __update_freq);


			const vector<OpenMM::Vec3>& minimized_positions = __system_topology.get_positions_in_nm();

			// check if minimization failed
			if (std::isnan(minimized_positions[0][0]))
				throw MinimizationError("die : minimization failed (in loop)");

			dbgmsg("minimized_positions = " << minimized_positions);
#ifndef NDEBUG
			const vector<OpenMM::Vec3>& forces = __system_topology.get_forces();
			dbgmsg("forces after minimization = " << forces);
#endif			
			// check if positions have converged
			double max_error = 0;
			for (int i = 0; i < initial_positions.size(); ++i) {
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
		cout << "Minimized in " << iter << " iterations, which took " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
	}

	void Modeler::init_openmm_positions() {
#ifndef NDEBUG
		dbgmsg("init_openmm_positions");
		for (auto &point : __positions) dbgmsg(point);
#endif
		__system_topology.init_positions(__positions);
	}

	void Modeler::init_openmm() {
		__system_topology.set_forcefield(*__ffield);
		__system_topology.init_particles(__topology);
		__system_topology.init_bonded(__topology, __use_constraints);
		__system_topology.init_integrator(__step_size_in_ps);

		if (__fftype == "kb") {

			// do nothing since force will be initialized when minimization starts
			
		} else if (__fftype == "phy") {
			
			__system_topology.init_physics_based_force(__topology);
			
		} else {
			throw Error("die : unsupported forcefield");
		}
		
	}

};
