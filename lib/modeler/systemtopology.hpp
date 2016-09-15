#ifndef SYSTEMTOPOLOGY_H
#define SYSTEMTOPOLOGY_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/geom3d.hpp"
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "pdbreader/molecule.hpp"
#include "modeler/topology.hpp"
#include "OpenMM.h"
using namespace std;

namespace Molib {
	class Atom;
	class Molecule;
};

namespace OMMIface {
	class ForceField;

	class SystemTopology {
	public:
		typedef enum {torsional=1, nonbond=2} options;
	private:
		OpenMM::System *system;
	    OpenMM::Integrator *integrator;
	    OpenMM::Context *context;

		OpenMM::HarmonicBondForce *bondStretch;
		OpenMM::HarmonicAngleForce *bondBend;
		OpenMM::PeriodicTorsionForce *bondTorsion;

		const ForceField *__ffield;
		
		int __kbforce_idx;
		vector<bool> masked;
		vector<double> masses;

		class AtomPoint {
		private:
			const Geom3D::Point __crd;
			Molib::Atom &__atom; 
		public:
			AtomPoint(const Geom3D::Point &crd, Molib::Atom &atom) : __crd(crd), __atom(atom) {}
			const Geom3D::Point& crd() const { return __crd; }
			Molib::Atom& get_atom() { return __atom; }
			void distance(double d) const {} // just dummy : needed by grid
			
			typedef vector<unique_ptr<AtomPoint>> UPVec;
			typedef vector<AtomPoint*> PVec;
			typedef ::Grid<AtomPoint> Grid;
		};
		
		struct ForceData {
			int force_idx, idx1, idx2, idx3, idx4;
			double length, angle;
			int periodicity;
			double phase, k;
		};
		vector<vector<ForceData>> bondStretchData, bondBendData, bondTorsionData;

	public:
		SystemTopology() : system(nullptr), integrator(nullptr), context(nullptr), __kbforce_idx(-1) {}
		~SystemTopology();
		static void loadPlugins();
		void mask(Topology &topology, const Molib::Atom::Vec &atoms);
		void unmask(Topology &topology, const Molib::Atom::Vec &atoms);

		void mask_forces(const int atom_idx, const set<int> &substruct);
		void unmask_forces(const int atom_idx, const set<int> &substruct);

		void init_integrator(const double step_size_in_ps);
		void init_particles(Topology &topology);
		void update_knowledge_based_force(Topology &topology, const vector<OpenMM::Vec3> &positions, const double dist_cutoff);
		void init_physics_based_force(Topology &topology);
		void init_bonded(Topology &topology, const bool use_constraints);
		void init_positions(const Geom3D::Point::Vec &crds);
		//~ const vector<OpenMM::Vec3>& get_positions_in_nm();
		vector<OpenMM::Vec3> get_positions_in_nm();
		vector<OpenMM::Vec3> get_forces();
		void minimize(const double tolerance, const double max_iterations);
		void set_forcefield(const ForceField &ffield) { __ffield = &ffield; }
	};
};
#endif
