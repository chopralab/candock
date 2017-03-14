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
#include "/usr/local/openmm/include/KBForce.h"
#include "topology.hpp"
#include "systemtopology.hpp"

using namespace std;

namespace Molib {
	class Score;
};

namespace OMMIface {
	
	class Modeler {
	public:
		class MinimizationError : public Error {
		public: 
			MinimizationError(const string &msg) : Error(msg) {}
		};
	private:
		const ForceField *__ffield;
		string __fftype;
		double __dist_cutoff_in_nm;
		double __tolerance;
		int __max_iterations;
		int __update_freq;
		double __position_tolerance_in_nm;
		bool __use_constraints;
		double __step_size_in_ps;

		Geom3D::Point::Vec __positions;
		Topology __topology;
		SystemTopology __system_topology;

	public:

		Modeler(const ForceField &ffield, const string &fftype, double dist_cutoff,
			double tolerance, int max_iterations, int update_freq, double position_tolerance,
			bool use_constraints, double step_size_in_fs) 
			: 
			__ffield(&ffield), __fftype(fftype), __dist_cutoff_in_nm(dist_cutoff * OpenMM::NmPerAngstrom),
			__tolerance(tolerance), __max_iterations(max_iterations), __update_freq(update_freq), 
			__position_tolerance_in_nm(position_tolerance * OpenMM::NmPerAngstrom), 
			__use_constraints(use_constraints), __step_size_in_ps(step_size_in_fs * OpenMM::PsPerFs)
			{}

		void mask(const Molib::Atom::Vec &atoms);
		void unmask(const Molib::Atom::Vec &atoms);
		
		void add_topology(const Molib::Atom::Vec &atoms);
		void add_crds(const Molib::Atom::Vec &atoms, const Geom3D::Point::Vec &crds);
		void add_random_crds(const Molib::Atom::Vec &atoms);

		Geom3D::Point::Vec get_state(const Molib::Atom::Vec &atoms);

		void minimize_state(const Molib::Molecule &ligand, const Molib::Molecule &receptor, const Molib::Score &score);
		void minimize_knowledge_based(const Molib::Molecule &ligand, const Molib::Molecule &receptor, const Molib::Score &score);
		void minimize_physical();

		void init_openmm_positions();
		void init_openmm();
		
		void set_max_iterations(const int max_iterations) { __max_iterations = max_iterations; }

	};
	
}

#endif
