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

using namespace std;

namespace OMMIface {
	
	class Modeler {
		ForceField &__ffield;
		string __fftype;
		double __dist_cutoff;
		bool __use_constraints;
		double __step_size_in_fs;
		
		map<void*, Topology> __topologies;
		map<void*, Geom3D::Point::Vec> __coords;

		class SystemTopology {
		public:
			typedef enum {torsional=1, non_bond=2} options;
		private:
			OpenMM::System *system;
		    OpenMM::Integrator*     integrator;
		    OpenMM::Context*  context;

			OpenMM::NonbondedForce nonbond;
			OpenMM::HarmonicBondForce bondStretch;
			OpenMM::HarmonicAngleForce bondBend;
			OpenMM::PeriodicTorsionForce bondTorsion;
			KBPlugin::KBForce *kbforce;
			
			vector<OpenMM::Vec3> initialPosInNm;
			map<const Molib::Atom*, const int> atom_to_index;

		public:
			void Modeler::SystemTopology::init_bonded(Topology &topology);
		};

	public:
		
		void add_topology(void *id, Molib::Atom::Vec &atoms);
		void add_coords(void *id, Geom3D::Point::Vec &coords);

		Geom3D::Point::Vec get_coords(void *id);
		Molib::Atom::Vec get_atoms(void *id);
		Molib::Atom::Vec get_atoms();

		void minimize_state();

		void set_forcefield(ForceField &ffield) { __ffield = ffield; }
		void set_forcefield_type(const string &fftype) { __fftype = fftype; }
		void set_distance_cutoff(const double dist_cutoff) { __dist_cutoff = dist_cutoff; }
		void set_use_constraints(const bool use_constraints) { __use_constraints = use_constraints; }
		void set_step_size_in_fs(const double step_size_in_fs) { __step_size_in_fs = step_size_in_fs; }

	};
	
}

#endif
