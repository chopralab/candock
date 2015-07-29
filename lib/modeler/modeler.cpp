#include "modeler.hpp"
#include "forcefield.hpp"
#include "moleculeinfo.hpp"
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

	void Modeler::loadPlugins() {
		// Load all available OpenMM plugins from their default location.
		
		dbgmsg("before loading plugins");
		OpenMM::Platform::loadPluginsFromDirectory
			(OpenMM::Platform::getDefaultPluginsDirectory());
		registerKBReferenceKernelFactories();
		dbgmsg("after loading plugins");
		
	}

	void Modeler::add_topology(void *id, Molib::Atom::Vec &atoms) {
		__topologies[id].init_topology(atoms, __ffield); // should also calculated bonded exclusions
	}

	void Modeler::add_coords(void *id, Geom3D::Point::Vec &coords) {
		__coords[id] = coords;
	}

	//~ void Modeler::SystemTopology::init_bonded(map<void*, Topology> &topologies) {
	void Modeler::SystemTopology::init_bonded(Topology &topology) {
			
		int warn = 0;
		
		dbgmsg("initializing openmm");
		
		// Specify the atoms and their properties:
		//  (1) System needs to know the masses.
		//  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
		//  (3) Collect default positions for initializing the simulation later.
		int idx = 0;
		for (auto &patom : topology.atoms) {
			const Molib::Atom &atom = *patom;
			if (!atom_to_index.count(&atom)) { // e.g., overlapping atoms from different segments
				atom_to_index.insert({&atom, idx++});
				const int type = topology.atom_to_type.at(&atom);
				try {
					const ForceField::AtomType& atype = __ffield.get_atom_type(type);
					system.addParticle(atype.mass);
					
					if (nonbond) nonbond->addParticle(atype.charge, atype.sigma, atype.epsilon);
	
					dbgmsg("add particle type = " << type << " crd = " << atom.crd() << " mass = " 
						<< atype.mass << " charge = " << atype.charge << " sigma = " << atype.sigma
						<< " epsilon = " << atype.epsilon << " representing atom = "
						<< atom << " at index = " << atom_to_index.at(&atom));
					const OpenMM::Vec3 posInNm(atom.crd().x() * OpenMM::NmPerAngstrom,
												atom.crd().y() * OpenMM::NmPerAngstrom,
												atom.crd().z() * OpenMM::NmPerAngstrom);
					initialPosInNm.push_back(posInNm);
				} catch (ParameterError& e) {
					cerr << e.what() << " (" << ++warn << ")" << endl;
				}
			}
		}
		// Process the bonds:
		//  (1) If we're using constraints, tell System about constrainable bonds;
		//      otherwise, tell HarmonicBondForce the bond stretch parameters 
		//      (tricky units!).
		//  (2) Create a list of bonds for generating nonbond exclusions.
		vector< pair<int,int> > bondPairs;
		for (auto &bond : topology.bonds) {
			dbgmsg("checkpoint0");
			const Molib::Atom &atom1 = *bond.first;
			dbgmsg(atom1);
			const Molib::Atom &atom2 = *bond.second;
			dbgmsg(atom2);
			const int idx1 = atom_to_index.at(&atom1);
			dbgmsg(idx1);
			const int idx2 = atom_to_index.at(&atom2);
			dbgmsg(idx2);
			const int type1 = topology.atom_to_type.at(&atom1);
			dbgmsg(type1);
			const int type2 = topology.atom_to_type.at(&atom2);
			dbgmsg(type2);
			try {				
				const ForceField::BondType& btype = __ffield.get_bond_type(type1, type2);
				if (__use_constraints && btype.can_constrain) { // Should we constrain C-H bonds?
					system.addConstraint(idx1, idx2, btype.length);
				} else {
					// Note factor of 2 for stiffness below because Amber specifies the constant
					// as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
					// it as used in the force term kx, with energy kx^2/2.
					if (bondStretch) bondStretch->addBond(idx1, idx2, btype.length, btype.k);
				}
				bondPairs.push_back({idx1, idx2});
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
		}
		dbgmsg("checkpoint3");
		// Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms.
		if (nonbond)
			nonbond->createExceptionsFromBonds(bondPairs, __ffield.coulomb14scale, __ffield.lj14scale);
		// Create the 1-2-3 bond angle harmonic terms.
		for (auto &angle : topology.angle) {
			dbgmsg("checkpoint4");
			const Molib::Atom &atom1 = *get<0>(angle);
			const Molib::Atom &atom2 = *get<1>(angle);
			const Molib::Atom &atom3 = *get<2>(angle);
			const int idx1 = atom_to_index.at(&atom1);
			const int idx2 = atom_to_index.at(&atom2);
			dbgmsg("checkpoint5");
			const int idx3 = atom_to_index.at(&atom3);
			const int type1 = topology.atom_to_type.at(&atom1);
			const int type2 = topology.atom_to_type.at(&atom2);
			const int type3 = topology.atom_to_type.at(&atom3);
			dbgmsg("checkpoint6");
			try {
				dbgmsg("determining angle type between atoms : " << endl
					<< atom1 << endl << atom2 << endl << atom3);
				const ForceField::AngleType& atype = 
					__ffield.get_angle_type(type1, type2, type3);
				if (bondBend) bondBend->addAngle(idx1, idx2, idx3, atype.angle, atype.k);
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
			dbgmsg("checkpoint8");
		}
		// Create the 1-2-3-4 bond torsion (dihedral) terms.
		for (auto &dihedral : topology.dihedral) {
			dbgmsg("checkpoint9");
			const Molib::Atom &atom1 = *get<0>(dihedral);
			const Molib::Atom &atom2 = *get<1>(dihedral);
			const Molib::Atom &atom3 = *get<2>(dihedral);
			const Molib::Atom &atom4 = *get<3>(dihedral);
			const int idx1 = atom_to_index.at(&atom1);
			const int idx2 = atom_to_index.at(&atom2);
			const int idx3 = atom_to_index.at(&atom3);
			const int idx4 = atom_to_index.at(&atom4);
			const int type1 = topology.atom_to_type.at(&atom1);
			const int type2 = topology.atom_to_type.at(&atom2);
			const int type3 = topology.atom_to_type.at(&atom3);
			const int type4 = topology.atom_to_type.at(&atom4);
			try {
				const ForceField::TorsionTypeVec& v_ttype = 
					__ffield.get_dihedral_type(type1, type2, type3, type4); // cannot make it const ??
				for (auto &ttype : v_ttype) {
					if (bondTorsion)
						bondTorsion->addTorsion(idx1, idx2, idx3, idx4, 
							ttype.periodicity, 
							ttype.phase,
							ttype.k);
				}
			} catch (ParameterError& e) {
				cerr << e.what() << " (" << ++warn << ")" << endl;
			}
		}
		// Create the 1-2-3-4 improper terms where 3 is the central atom
		for (auto &dihedral : topology.improper) {
			dbgmsg("checkpoint9");
			const Molib::Atom &atom1 = *get<0>(dihedral);
			const Molib::Atom &atom2 = *get<1>(dihedral);
			const Molib::Atom &atom3 = *get<2>(dihedral);
			const Molib::Atom &atom4 = *get<3>(dihedral);
			const int idx1 = atom_to_index.at(&atom1);
			const int idx2 = atom_to_index.at(&atom2);
			const int idx3 = atom_to_index.at(&atom3);
			const int idx4 = atom_to_index.at(&atom4);
			const int type1 = topology.atom_to_type.at(&atom1);
			const int type2 = topology.atom_to_type.at(&atom2);
			const int type3 = topology.atom_to_type.at(&atom3);
			const int type4 = topology.atom_to_type.at(&atom4);
			try {
				const ForceField::TorsionTypeVec& v_ttype = 
					__ffield.get_improper_type(type1, type2, type3, type4); // cannot make it const ??
				for (auto &ttype : v_ttype) {
					if (bondTorsion)
						bondTorsion->addTorsion(idx1, idx2, idx3, idx4, 
							ttype.periodicity, 
							ttype.phase,
							ttype.k);
				}
			} catch (ParameterError& e) {
				dbgmsg(e.what() << " (WARNINGS ARE NOT INCREASED)");
			}
		}
		// Create the knowledge-based forcefield terms.
		if (kbforce) {
			kbforce->setStep(__ffield.step);
			kbforce->setScale(__ffield.scale);
			dbgmsg("kbforce scale = " << __ffield.scale 
				<< " step = " << __ffield.step);
			for (auto &bond : topology.kbforce) {
				const Molib::Atom &atom1 = *bond.first;
				const Molib::Atom &atom2 = *bond.second;
				const int idx1 = atom_to_index.at(&atom1);
				const int idx2 = atom_to_index.at(&atom2);
				const int type1 = topology.atom_to_type.at(&atom1);
				const int type2 = topology.atom_to_type.at(&atom2);
				try {
					const ForceField::KBType& kbtype = 
						__ffield.get_kb_force_type(atom1, atom2, type1, type2);
					kbforce->addBond(idx1, idx2, 
						const_cast<vector<double>&>(kbtype.potential), 
						const_cast<vector<double>&>(kbtype.derivative));
				} catch (ParameterError& e) {
					cerr << e.what() << " (" << ++warn << ")" << endl;
				}
			}
		}
		// Choose an Integrator for advancing time, and a Context connecting the
		// System with the Integrator for simulation. Let the Context choose the
		// best available Platform. Initialize the configuration from the default
		// positions we collected above. Initial velocities will be zero.
		omm->integrator = new OpenMM::VerletIntegrator(__step_size_in_fs * OpenMM::PsPerFs);
		omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
		omm->context->setPositions(initialPosInNm);
		dbgmsg("REMARK  Using OpenMM platform "	<< omm->context->getPlatform().getName());
		if (warn > 0) {
			throw Error("die : missing parameters detected for molecule "
				+ __ligand.name());
		}
	}
	pair<Molib::Molecule, Molib::Molecule> OMM::__get_openmm_state(const MyOpenMMData* omm, 
		const Molib::Molecule &receptor, const Molib::Molecule &ligand) const {
		int infoMask = OpenMM::State::Positions;
		// Forces are also available (and cheap).
		const OpenMM::State state = omm->context->getState(infoMask);
		// Copy OpenMM positions into atoms array and change units from nm to Angstroms.
		const vector<OpenMM::Vec3>& positions_in_nm = state.getPositions();
		Molib::Molecule minimized_receptor(receptor); // receptor and ligand are two separate molecules
		Molib::Molecule minimized_ligand(ligand); // receptor and ligand are two separate molecules

		int i = 0;
		for (auto &patom : minimized_receptor.get_atoms()) {
			 // change coordinates of initial system to new coordinates
			patom->set_crd(Geom3D::Coordinate(positions_in_nm[i][0] * OpenMM::AngstromsPerNm,
										positions_in_nm[i][1] * OpenMM::AngstromsPerNm,
										positions_in_nm[i][2] * OpenMM::AngstromsPerNm));
			++i;
		}
		for (auto &patom : minimized_ligand.get_atoms()) {
			 // change coordinates of initial system to new coordinates
			patom->set_crd(Geom3D::Coordinate(positions_in_nm[i][0] * OpenMM::AngstromsPerNm,
										positions_in_nm[i][1] * OpenMM::AngstromsPerNm,
										positions_in_nm[i][2] * OpenMM::AngstromsPerNm));
			++i;
		}
		return {minimized_receptor, minimized_ligand};
	}
	double OMM::__get_openmm_energy(const MyOpenMMData* omm) const {
		/* Returns energy (and time)
		 * 
		 */
		int infoMask = 0;
		
		//~ infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
		infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
		const OpenMM::State state = omm->context->getState(infoMask);
		double time_in_ps = state.getTime(); // OpenMM time is in ps already
		double energy_in_kcal = (state.getPotentialEnergy() + state.getKineticEnergy())
						* OpenMM::KcalPerKJ;
		
		return energy_in_kcal;
	}

	void OMM::minimize(double tolerance, int max_iterations, int update_freq) {
		Benchmark::reset();
		cout << "Doing energy minimization of ligand " << __ligand.name() << endl;

		if (__fftype == "phy") {
			OpenMM::LocalEnergyMinimizer::minimize(*__omm->context, 
								tolerance, max_iterations);
		} else {
			// for knowledge-based forcefield we implement a custom update nonbond
			// function
			
			int iter = 0;
	
			// get initial positions
			OpenMM::State initial_state = __omm->context->getState(OpenMM::State::Positions);
			vector<OpenMM::Vec3> initial_positions = initial_state.getPositions();
	
			while (iter < max_iterations) {
				
				dbgmsg("starting minimization step = " << iter);
				
				OpenMM::LocalEnergyMinimizer::minimize(*__omm->context, 
											tolerance, update_freq);
	
				/**
				 *  update nonbond list
				 */
				 
				 // remove old knowledge-based force
				__omm->system->removeForce(__kbforce_idx);
				
				// recalculate distances
				auto ret = this->get_state(__receptor, __ligand);
				Molib::Molecule minimized_receptor(ret.first);
				Molib::Molecule minimized_ligand(ret.second);
				
				Topology topology;
				topology.get_topology(minimized_receptor, __ffield)
					.get_topology(minimized_ligand, __ffield)
					.get_kb_force_info(minimized_receptor, minimized_ligand, __dist_cutoff);				
	
				// prepare atom_to_index
				map<const Molib::Atom*, const int> atom_to_index;
				int idx = 0;
				for (auto &patom : topology.atom) {
					const Molib::Atom &atom = *patom;
					atom_to_index.insert({&atom, idx++});
				}
	
				// create new knowledge-based force
				KBPlugin::KBForce* kbforce = new KBPlugin::KBForce();
	
				kbforce->setStep(__ffield.step);
				kbforce->setScale(__ffield.scale);
				dbgmsg("kbforce scale = " << __ffield.scale 
					<< " step = " << __ffield.step);
				
				for (auto &bond : topology.kbforce) {
					const Molib::Atom &atom1 = *bond.first;
					const Molib::Atom &atom2 = *bond.second;
					const int idx1 = atom_to_index.at(&atom1);
					const int idx2 = atom_to_index.at(&atom2);
					const int type1 = topology.atom_to_type.at(&atom1);
					const int type2 = topology.atom_to_type.at(&atom2);
					try {
						const ForceField::KBType& kbtype = 
							__ffield.get_kb_force_type(atom1, atom2, type1, type2);
						kbforce->addBond(idx1, idx2, 
							const_cast<vector<double>&>(kbtype.potential), 
							const_cast<vector<double>&>(kbtype.derivative));
					} catch (ParameterError& e) {
						cerr << e.what() << endl;
					}
				}
				//~ dbgmsg("before adding new force numBonds = " << ((KBPlugin::KBForce&)__omm->system->getForce(__kbforce_idx)).getNumBonds()
					//~ << " __kbforce_idx = " << __kbforce_idx);
				dbgmsg("before adding new force");
				
				// and add it to the system
				__kbforce_idx = __omm->system->addForce(kbforce);
	
				dbgmsg("after adding new force numBonds = " << ((KBPlugin::KBForce&)__omm->system->getForce(__kbforce_idx)).getNumBonds()
					<< " __kbforce_idx = " << __kbforce_idx);
				// get positions after minimization
				OpenMM::State state = __omm->context->getState(OpenMM::State::Positions);
				vector<OpenMM::Vec3> minimized_positions = state.getPositions();
	
				dbgmsg("before checking if positions converged");
	
				// check if positions don't change too much
				double max_error = 0;
				for (int i = 0; i < __omm->system->getNumParticles(); ++i) {
					OpenMM::Vec3 dif_pos = initial_positions[i] - minimized_positions[i];
					 const double error = dif_pos.dot(dif_pos);
					 if (error > max_error)
						max_error = error;
					 
				}
	
				dbgmsg("after checking if positions converged");
				
				// convergence reached
				//~ if (sqrt(max_error) < pos_tol)
				if (sqrt(max_error) < 0.0001)
					break;
					
				iter += update_freq;						
	
				initial_positions = minimized_positions;
				dbgmsg("ending minimization step iter = " << iter);
	
			}
		}
		cout << "time to minimize ligand " << __ligand.name() << " took " 
			<< Benchmark::seconds_from_start() << " wallclock seconds" << endl;
	}

	void OMM::md(const double simulation_time_in_ps, 
					const double report_interval_in_fs,
					const bool want_energy) const {
		// Run the simulation:
		//  (1) Write the first line of the PDB file and the initial configuration.
		//  (2) Run silently entirely within OpenMM between reporting intervals.
		//  (3) Write a PDB frame when the time comes.
		const int num_silent_steps = (int)(report_interval_in_fs / __step_size_in_fs + 0.5);
		for (int frame=1; ; ++frame) {
			double time, energy;
			//~ Molib::Molecules mols = __get_openmm_state(__omm, want_energy, __receptor, time, energy); // MAKE IT WORK LATER!!!
			//~ __myWritePDBFrame(frame, time, energy, atoms);			
			if (time >= simulation_time_in_ps)
			break;
			__step_with_openmm(__omm, num_silent_steps);
		}
    }
	Components OMM::get_energy_components(const Molib::Molecule &receptor, 
		const Molib::Molecule &ligand, const double cur_dist_cutoff) {
		/* Get energy components : interaction energy, internal ligand energy, 
		 * ligand torsional energy, internal receptor energy, receptor torsional energy
		 * 
		 */
		assert(cur_dist_cutoff <= __dist_cutoff);
		
		Topology complex_info;
		complex_info.get_topology(receptor, __ffield)
			.get_topology(ligand, __ffield)
			.get_kb_force_info(receptor, ligand, cur_dist_cutoff);
		Topology interaction_info;
		interaction_info.get_topology(receptor, __ffield)
			.get_topology(ligand, __ffield)
			.get_interaction_force_info(receptor, ligand, cur_dist_cutoff);
		Topology ligand_info;					
		ligand_info.get_topology(ligand, __ffield)
			.get_kb_force_info(Molib::Molecule("dummy"), ligand, cur_dist_cutoff);
		Topology receptor_info;					
		receptor_info.get_topology(receptor, __ffield)
			.get_kb_force_info(receptor, Molib::Molecule("dummy"), cur_dist_cutoff);
		
		dbgmsg("after moleculeinfo");
		MyOpenMMData omm_tot, omm_ie, omm_lie, omm_lte, omm_rie, omm_rte;
		__initialize_openmm(&omm_tot, complex_info, OMMIface::non_bond|OMMIface::torsional); 
		double total_energy = __get_openmm_energy(&omm_tot);
		
		__initialize_openmm(&omm_ie, interaction_info, OMMIface::non_bond); 
		double interaction_energy = __get_openmm_energy(&omm_ie);
		
		__initialize_openmm(&omm_lie, ligand_info, OMMIface::non_bond); 
		double ligand_internal_energy = __get_openmm_energy(&omm_lie);
		
		__initialize_openmm(&omm_lte, ligand_info, OMMIface::torsional); 
		double ligand_torsional_energy = __get_openmm_energy(&omm_lte);
		
		__initialize_openmm(&omm_rie, receptor_info, OMMIface::non_bond); 
		double receptor_internal_energy = __get_openmm_energy(&omm_rie);
		
		__initialize_openmm(&omm_rte, receptor_info, OMMIface::torsional); 
		double receptor_torsional_energy = __get_openmm_energy(&omm_rte);
		
		 return Components(total_energy, interaction_energy, ligand_internal_energy, 
			ligand_torsional_energy, receptor_internal_energy, 
			receptor_torsional_energy);
	}
};
