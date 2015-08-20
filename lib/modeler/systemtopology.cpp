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
	extern "C" OPENMM_EXPORT void registerKBReferenceKernelFactories();

    SystemTopology::~SystemTopology() {
		dbgmsg("calling destructor of SystemTopology");
		delete context; delete integrator; delete system;
	}

	void SystemTopology::loadPlugins() {
		// Load all available OpenMM plugins from their default location.
		
		dbgmsg("before loading plugins");
		OpenMM::Platform::loadPluginsFromDirectory
			(OpenMM::Platform::getDefaultPluginsDirectory());
		registerKBReferenceKernelFactories();
		dbgmsg("after loading plugins");
		
	}

	void SystemTopology::mask(Topology &topology, const Molib::Atom::Vec &atoms) {
		for (auto &patom : atoms) {
			int idx = topology.get_index(*patom);
			system->setParticleMass(idx, 0);
			masked[idx] = true;
		}
	}
	
	void SystemTopology::unmask(Topology &topology, const Molib::Atom::Vec &atoms) {
		for (auto &patom : atoms) {
			int idx = topology.get_index(*patom);
			system->setParticleMass(idx, masses[idx]);
			masked[idx] = false;
		}
	}
	
	void SystemTopology::init_integrator(const double step_size_in_ps) {
		// Choose an Integrator for advancing time, and a Context connecting the
		// System with the Integrator for simulation. Let the Context choose the
		// best available Platform. Initialize the configuration from the default
		// positions we collected above. Initial velocities will be zero.
		integrator = new OpenMM::VerletIntegrator(step_size_in_ps);
		context = new OpenMM::Context(*system, *integrator);
		dbgmsg("REMARK  Using OpenMM platform "	<< context->getPlatform().getName());
	}

	void SystemTopology::init_particles(Topology &topology) {

		int warn = 0;

		system = new OpenMM::System();

		dbgmsg("topology.atoms.size() = " << topology.atoms.size());
		
		// Specify the atoms and their properties:
		//  (1) System needs to know the masses.
		//  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
		//  (3) Collect default positions for initializing the simulation later.
		for (auto &patom : topology.atoms) {
			const Molib::Atom &atom = *patom;
			const int type = topology.get_type(atom);
			try {
				const ForceField::AtomType& atype = __ffield->get_atom_type(type);
				system->addParticle(atype.mass);

				masses.push_back(atype.mass);
				masked.push_back(true); // at the start mask everything
				
				dbgmsg("add particle type = " << type << " crd = " << atom.crd() << " mass = " 
					<< atype.mass << " charge = " << atype.charge << " sigma = " << atype.sigma
					<< " epsilon = " << atype.epsilon << " representing atom = "
					<< atom << " at index = " << topology.get_index(atom));
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
		}
		if (warn > 0) {
			throw Error("die : missing parameters detected");
		}
	}

	void SystemTopology::update_knowledge_based_force(Topology &topology, const vector<OpenMM::Vec3> &positions, const double dist_cutoff) {

		// we have to set new coordinates for atoms
		AtomPoint::UPVec atompoints;
		for (int i = 0; i < positions.size(); ++i) {
			atompoints.push_back(unique_ptr<AtomPoint>(new AtomPoint(
				Geom3D::Point(positions[i][0], positions[i][1], positions[i][2]),
				*topology.atoms[i])));
		}

		dbgmsg("atompoints.size() = " << atompoints.size());
		
		// generate grid
		AtomPoint::Grid grid(atompoints);
		
		// knowledge-based forces between atoms except between bonded exclusions
		Molib::Atom::Set visited;
		Topology::Bonds nonbondlist;

		for (auto &pa1 : atompoints) {
			Molib::Atom &atom1 = pa1->get_atom();
			visited.insert(&atom1);
			for (auto &pa2 : grid.get_neighbors(pa1->crd(), dist_cutoff)) {
				Molib::Atom &atom2 = pa2->get_atom();
				if (!visited.count(&atom2)) {
					if (!topology.bonded_exclusions.count({&atom1, &atom2})) {
						nonbondlist.push_back({&atom1, &atom2});
					}
				}
			}
		}
		
		dbgmsg("NONBOND LIST (KBFORCES) : " << endl << nonbondlist);

		// remove force if it exists
		if (__kbforce_idx != -1)
			system->removeForce(__kbforce_idx);

		KBPlugin::KBForce* kbforce = new KBPlugin::KBForce();
		kbforce->setStep(__ffield->step);
		kbforce->setScale(__ffield->scale);
		dbgmsg("kbforce scale = " << __ffield->scale << " step = " << __ffield->step);

		int warn = 0;

		// init OpenMM's kbforce object
		for (auto &bond : nonbondlist) {
			const Molib::Atom &atom1 = *bond.first;
			const Molib::Atom &atom2 = *bond.second;
			const int idx1 = topology.get_index(atom1);
			const int idx2 = topology.get_index(atom2);
			const int type1 = topology.get_type(atom1);
			const int type2 = topology.get_type(atom2);
			if (!masked[idx1] && !masked[idx2]) { // don't make the force if one or both atoms are masked
				try {
					const ForceField::KBType& kbtype = 
						__ffield->get_kb_force_type(atom1, atom2, type1, type2);
					kbforce->addBond(idx1, idx2, 
						const_cast<vector<double>&>(kbtype.potential), 
						const_cast<vector<double>&>(kbtype.derivative));
					dbgmsg("adding kbforce between idx1 = " << idx1 << " and idx2 = " << idx2);
				} catch (ParameterError& e) {
					cerr << e.what() << " (" << ++warn << ")" << endl;
				}
			}
		}	
		// and add it back to the system
		__kbforce_idx = system->addForce(kbforce);
		context->reinitialize();
		context->setPositions(positions);
		if (warn > 0) {
			throw Error("die : missing parameters detected");
		}

	}
	
	void SystemTopology::init_physics_based_force(Topology &topology) {

		int warn = 0;

		OpenMM::NonbondedForce *nonbond = new OpenMM::NonbondedForce();
 		system->addForce(nonbond);

		for (auto &patom : topology.atoms) {
			const Molib::Atom &atom = *patom;
			const int type = topology.get_type(atom);
			try {
				const ForceField::AtomType& atype = __ffield->get_atom_type(type);

				nonbond->addParticle(atype.charge, atype.sigma, atype.epsilon);

				dbgmsg("add particle type = " << type << " crd = " << atom.crd() << " mass = " 
					<< atype.mass << " charge = " << atype.charge << " sigma = " << atype.sigma
					<< " epsilon = " << atype.epsilon << " representing atom = "
					<< atom << " at index = " << topology.get_index(atom));

			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
		}
		
		vector< pair<int,int> > bondPairs;
		for (auto &bond : topology.bonds) {
			dbgmsg("checkpoint0");
			const Molib::Atom &atom1 = *bond.first;
			dbgmsg(atom1);
			const Molib::Atom &atom2 = *bond.second;
			dbgmsg(atom2);
			const int idx1 = topology.get_index(atom1);
			dbgmsg(idx1);
			const int idx2 = topology.get_index(atom2);
			dbgmsg(idx2);
			bondPairs.push_back({idx1, idx2});
		}
		dbgmsg("checkpoint3");
		// Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms.
		nonbond->createExceptionsFromBonds(bondPairs, __ffield->coulomb14scale, __ffield->lj14scale);

		if (warn > 0) {
			throw Error("die : missing parameters detected");
		}
		
	}
	
	void SystemTopology::init_bonded(Topology &topology, const bool use_constraints) {
			
		int warn = 0;

		OpenMM::HarmonicBondForce *bondStretch = new OpenMM::HarmonicBondForce();
		OpenMM::HarmonicAngleForce *bondBend = new OpenMM::HarmonicAngleForce();
		OpenMM::PeriodicTorsionForce *bondTorsion = new OpenMM::PeriodicTorsionForce();

		system->addForce(bondStretch);
		system->addForce(bondBend);
		system->addForce(bondTorsion);

		dbgmsg("initializing openmm");
		
		// Process the bonds:
		//  (1) If we're using constraints, tell System about constrainable bonds;
		//      otherwise, tell HarmonicBondForce the bond stretch parameters 
		//      (tricky units!).
		//  (2) Create a list of bonds for generating nonbond exclusions.
		for (auto &bond : topology.bonds) {
			dbgmsg("checkpoint0");
			const Molib::Atom &atom1 = *bond.first;
			dbgmsg(atom1);
			const Molib::Atom &atom2 = *bond.second;
			dbgmsg(atom2);
			const int idx1 = topology.get_index(atom1);
			dbgmsg(idx1);
			const int idx2 = topology.get_index(atom2);
			dbgmsg(idx2);
			const int type1 = topology.get_type(atom1);
			dbgmsg(type1);
			const int type2 = topology.get_type(atom2);
			dbgmsg(type2);
			try {				
				const ForceField::BondType& btype = __ffield->get_bond_type(type1, type2);
				if (use_constraints && btype.can_constrain) { // Should we constrain C-H bonds?
					system->addConstraint(idx1, idx2, btype.length);
				} else {
					// Note factor of 2 for stiffness below because Amber specifies the constant
					// as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
					// it as used in the force term kx, with energy kx^2/2.
					bondStretch->addBond(idx1, idx2, btype.length, btype.k);
				}
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
		}
		dbgmsg("checkpoint3");

		// Create the 1-2-3 bond angle harmonic terms.
		for (auto &angle : topology.angles) {
			dbgmsg("checkpoint4");
			const Molib::Atom &atom1 = *get<0>(angle);
			const Molib::Atom &atom2 = *get<1>(angle);
			const Molib::Atom &atom3 = *get<2>(angle);
			const int idx1 = topology.get_index(atom1);
			const int idx2 = topology.get_index(atom2);
			dbgmsg("checkpoint5");
			const int idx3 = topology.get_index(atom3);
			const int type1 = topology.get_type(atom1);
			const int type2 = topology.get_type(atom2);
			const int type3 = topology.get_type(atom3);
			dbgmsg("checkpoint6");
			try {
				dbgmsg("determining angle type between atoms : " << endl
					<< atom1 << endl << atom2 << endl << atom3);
				const ForceField::AngleType& atype = 
					__ffield->get_angle_type(type1, type2, type3);
				bondBend->addAngle(idx1, idx2, idx3, atype.angle, atype.k);
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
			dbgmsg("checkpoint8");
		}
		
		// Create the 1-2-3-4 bond torsion (dihedral) terms.
		for (auto &dihedral : topology.dihedrals) {
			dbgmsg("checkpoint9");
			const Molib::Atom &atom1 = *get<0>(dihedral);
			const Molib::Atom &atom2 = *get<1>(dihedral);
			const Molib::Atom &atom3 = *get<2>(dihedral);
			const Molib::Atom &atom4 = *get<3>(dihedral);
			const int idx1 = topology.get_index(atom1);
			const int idx2 = topology.get_index(atom2);
			const int idx3 = topology.get_index(atom3);
			const int idx4 = topology.get_index(atom4);
			const int type1 = topology.get_type(atom1);
			const int type2 = topology.get_type(atom2);
			const int type3 = topology.get_type(atom3);
			const int type4 = topology.get_type(atom4);
			try {
				const ForceField::TorsionTypeVec& v_ttype = 
					__ffield->get_dihedral_type(type1, type2, type3, type4); // cannot make it const ??
				for (auto &ttype : v_ttype) {
					bondTorsion->addTorsion(idx1, idx2, idx3, idx4, ttype.periodicity, ttype.phase,	ttype.k);
				}
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
		}
		
		// Create the 1-2-3-4 improper terms where 3 is the central atom
		for (auto &dihedral : topology.impropers) {
			dbgmsg("checkpoint9");
			const Molib::Atom &atom1 = *get<0>(dihedral);
			const Molib::Atom &atom2 = *get<1>(dihedral);
			const Molib::Atom &atom3 = *get<2>(dihedral);
			const Molib::Atom &atom4 = *get<3>(dihedral);
			const int idx1 = topology.get_index(atom1);
			const int idx2 = topology.get_index(atom2);
			const int idx3 = topology.get_index(atom3);
			const int idx4 = topology.get_index(atom4);
			const int type1 = topology.get_type(atom1);
			const int type2 = topology.get_type(atom2);
			const int type3 = topology.get_type(atom3);
			const int type4 = topology.get_type(atom4);
			try {
				const ForceField::TorsionTypeVec& v_ttype = 
					__ffield->get_improper_type(type1, type2, type3, type4); // cannot make it const ??
				for (auto &ttype : v_ttype) {
					bondTorsion->addTorsion(idx1, idx2, idx3, idx4, ttype.periodicity, ttype.phase, ttype.k);
				}
			} catch (ParameterError& e) {
				dbgmsg(e.what() << " (WARNINGS ARE NOT INCREASED)");
			}
		}
		if (warn > 0) {
			throw Error("die : missing parameters detected");
		}
	}

	void SystemTopology::init_positions(const Geom3D::Point::Vec &crds) { 
		dbgmsg("entering init_positions crds.size() = " << crds.size());
		vector<OpenMM::Vec3> positions_in_nm;
		for (auto &crd: crds) {
			positions_in_nm.push_back(OpenMM::Vec3(crd.x() * OpenMM::NmPerAngstrom, 
				crd.y() * OpenMM::NmPerAngstrom, crd.z() * OpenMM::NmPerAngstrom));
		}
		
		context->setPositions(positions_in_nm);
		
		dbgmsg("exiting init_positions");
	}

	vector<OpenMM::Vec3> SystemTopology::get_positions_in_nm() {
		return context->getState(OpenMM::State::Positions).getPositions();		
	}

	void SystemTopology::minimize(const double tolerance, const double max_iterations) {
		OpenMM::LocalEnergyMinimizer::minimize(*context, tolerance, max_iterations);
	}

};
