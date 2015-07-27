#ifndef OMM_H
#define OMM_H
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
using namespace std;
/* 
 * this is a layer between Molib and OpenMM
 */

namespace OMMIface {
	class ForceField;
	class MoleculeInfo;
	typedef enum {torsional=1, non_bond=2} energy_options;

	typedef tuple<double, double, double, double, double, double> Components;
	typedef map<Molib::Molecule*, Components> Energies;
	const string print_energies_title();
	class OMM {
		struct MyOpenMMData {
		    MyOpenMMData() : system(0), context(0), integrator(0) {}
		    ~MyOpenMMData() {
				dbgmsg("calling destructor of MyOpenMMData");
				delete context; delete integrator; delete system;
			}
		    OpenMM::System*         system;
		    OpenMM::Integrator*     integrator;
		    OpenMM::Context*  context;
		};
		void __initialize_openmm(MyOpenMMData*, const MoleculeInfo&, unsigned int);
		pair<Molib::Molecule, Molib::Molecule> __get_openmm_state(const MyOpenMMData* omm, 
			const Molib::Molecule &receptor, const Molib::Molecule &ligand) const;
		double __get_openmm_energy(const MyOpenMMData* omm) const;
		void __step_with_openmm(const MyOpenMMData* omm, int numSteps) const { omm->integrator->step(numSteps);	}
		//~ void __terminate_openmm(MyOpenMMData* omm) { delete omm; }
		MyOpenMMData *__omm;
		const ForceField &__ffield;
		const string &__fftype;
		const double __dist_cutoff;
		const bool __use_constraints;
		const double __step_size_in_fs; // integration step size (fs)
		const Molib::Molecule &__ligand;
	public:
		OMM(const Molib::Molecule &receptor, const Molib::Molecule &ligand,
			const ForceField &ff, const string &fftype, const double dist_cutoff, 
			const bool use_constraints=false, const double step_size_in_fs=2);
		~OMM() {
			//~ dbgmsg("calling destructor of OMM");
			//~ cout << "calling destructor of OMM for molecule " << __ligand.name() << endl;
			//~ __terminate_openmm(__omm); // Clean up OpenMM data structures.
			delete __omm; // Clean up OpenMM data structures.
		}
		void minimize(double tolerance, int max_iterations, int update_freq);
		//~ void minimize(double tolerance, int max_iterations, int update_freq);
		void md(const double=100, const double=10, const bool=false) const;
		pair<Molib::Molecule, Molib::Molecule> get_state(
			const Molib::Molecule &receptor, const Molib::Molecule &ligand) const { 
			return __get_openmm_state(__omm, receptor, ligand);
		}
	static void loadPlugins();
	Components get_energy_components(const Molib::Molecule &receptor, 
		const Molib::Molecule &ligand, const double cur_dist_cutoff);
	};
};

ostream& operator<<(ostream& os, const OMMIface::Energies& energies);

#endif
