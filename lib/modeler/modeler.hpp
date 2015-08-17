#ifndef MODELER_H
#define MODELER_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/coordinate.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"
#include "helper/benchmark.hpp"
#include "OpenMM.h"
#include "kbforce/openmmapi/include/KBForce.h"
#include "topology.hpp"
#include "systemtopology.hpp"

using namespace std;

namespace OMMIface {
	
	class Modeler {
		ForceField *__ffield;
		string __fftype;
		double __dist_cutoff_in_nm;
		bool __use_constraints;
		double __step_size_in_fs;
		double __tolerance;
		double __position_tolerance_in_nm;
		int __max_iterations;
		int __update_freq;

		Geom3D::Point::Vec __positions;
		Topology __topology;
		SystemTopology __system_topology;

	public:

		void mask(const Molib::Atom::Vec &atoms);
		void unmask(const Molib::Atom::Vec &atoms);
		
		void add_topology(const Molib::Atom::Vec &atoms);
		void add_crds(const Molib::Atom::Vec &atoms, const Geom3D::Point::Vec &crds);

		Geom3D::Point::Vec get_state(const Molib::Atom::Vec &atoms);

#ifndef NDEBUG
		void minimize_state(Molib::Molecule &ligand, Molib::Molecule &receptor);
		void minimize_knowledge_based(Molib::Molecule &ligand, Molib::Molecule &receptor);
#else
		void minimize_state() {
		void minimize_knowledge_based();
#endif
		void minimize_physical();

		void init_openmm_positions();
		void init_openmm();

		void set_forcefield(ForceField &ffield) { __ffield = &ffield; }
		void set_forcefield_type(const string &fftype) { __fftype = fftype; }
		void set_distance_cutoff(const double dist_cutoff) { __dist_cutoff_in_nm = dist_cutoff * OpenMM::NmPerAngstrom; }
		void set_use_constraints(const bool use_constraints) { __use_constraints = use_constraints; }
		void set_step_size_in_fs(const double step_size_in_fs) { __step_size_in_fs = step_size_in_fs; }
		void set_tolerance(const double tolerance) { __tolerance = tolerance; }
		void set_max_iterations(const int max_iterations) { __max_iterations = max_iterations; }
		void set_update_freq(const int update_freq) { __update_freq = update_freq; }
		void set_position_tolerance(const double position_tolerance) { __position_tolerance_in_nm = position_tolerance * OpenMM::NmPerAngstrom; }
	};
	
}

#endif
