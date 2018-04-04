#include "systemtopology.hpp"
#include "modeler.hpp"

#include <openmm/Platform.h>
#include <openmm/HarmonicBondForce.h>
#include <openmm/HarmonicAngleForce.h>
#include <openmm/PeriodicTorsionForce.h>
#include <openmm/VerletIntegrator.h>
#include <openmm/LangevinIntegrator.h>
#include <openmm/BrownianIntegrator.h>
#include <openmm/AndersenThermostat.h>
#include <openmm/NonbondedForce.h>
#include <openmm/CustomNonbondedForce.h>
#include <openmm/LocalEnergyMinimizer.h>

#include "forcefield.hpp"
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

namespace OMMIface
{

SystemTopology::~SystemTopology()
{
        dbgmsg("calling destructor of SystemTopology");
        delete context;
        delete integrator;
        delete system;
}

void SystemTopology::loadPlugins(const std::string &extra_dir)
{
        try
        {
                // Load all available OpenMM plugins from their default location.

                dbgmsg("before loading plugins");
                OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());

                if (!extra_dir.empty())
                {
                        OpenMM::Platform::loadPluginsFromDirectory(extra_dir);
                }

                dbgmsg("after loading plugins");
        }
        catch (const std::exception &e)
        {
                log_error << "The crash is " << e.what() << endl;
        }
}

void SystemTopology::mask(Topology &topology, const Molib::Atom::Vec &atoms)
{
        set<int> substruct;

        for (auto &patom : atoms)
        {
                int idx = topology.get_index(*patom);
                dbgmsg("masking particle idx = " << idx);
                substruct.insert(idx);
        }

        for (auto &idx : substruct)
        {
                system->setParticleMass(idx, 0);
                masked[idx] = true;

                mask_forces(idx, substruct);
        }

        bondStretch->updateParametersInContext(*context);
        bondBend->updateParametersInContext(*context);
        bondTorsion->updateParametersInContext(*context);
}

void SystemTopology::unmask(Topology &topology, const Molib::Atom::Vec &atoms)
{
        set<int> substruct;

        for (auto &patom : atoms)
        {
                int idx = topology.get_index(*patom);
                substruct.insert(idx);
        }

        for (auto &idx : substruct)
        {
                dbgmsg("unmasking particle idx = " << idx << " mass = " << masses[idx]);
                system->setParticleMass(idx, masses[idx]);
                masked[idx] = false;

                unmask_forces(idx, substruct);
        }

        bondStretch->updateParametersInContext(*context);
        bondBend->updateParametersInContext(*context);
        bondTorsion->updateParametersInContext(*context);
}

void SystemTopology::mask_forces(const int atom_idx, const set<int> &substruct)
{
        for (auto &data : bondStretchData[atom_idx])
        { // get all forces involving this atom's idx
                if (substruct.count(data.idx1) && substruct.count(data.idx2))
                        bondStretch->setBondParameters(data.force_idx, data.idx1, data.idx2, data.length, 0.0);
        }

        for (auto &data : bondBendData[atom_idx])
        { // get all forces involving this atom's idx
                if (substruct.count(data.idx1) && substruct.count(data.idx2) && substruct.count(data.idx3))
                        bondBend->setAngleParameters(data.force_idx, data.idx1, data.idx2, data.idx3, data.angle, 0.0);
        }

        for (auto &data : bondTorsionData[atom_idx])
        { // get all forces involving this atom's idx
                if (substruct.count(data.idx1) && substruct.count(data.idx2) && substruct.count(data.idx3) && substruct.count(data.idx4))
                {

                        int idx1, idx2, idx3, idx4;
                        int periodicity;
                        double phase, k;

                        bondTorsion->getTorsionParameters(data.force_idx, idx1, idx2, idx3, idx4, periodicity, phase, k);
                        dbgmsg("particles of force " << data.force_idx);
                        dbgmsg(idx1 << " " << data.idx1);
                        dbgmsg(idx2 << " " << data.idx2);
                        dbgmsg(idx3 << " " << data.idx3);
                        dbgmsg(idx4 << " " << data.idx4);
                        dbgmsg(periodicity << " " << data.periodicity);

                        if (idx1 != data.idx1)
                                throw Error("die : particles changed");

                        if (idx2 != data.idx2)
                                throw Error("die : particles changed");

                        if (idx3 != data.idx3)
                                throw Error("die : particles changed");

                        if (idx4 != data.idx4)
                                throw Error("die : particles changed");

                        bondTorsion->setTorsionParameters(data.force_idx, data.idx1, data.idx2, data.idx3, data.idx4, data.periodicity, data.phase, 0.0);
                }
        }
}

void SystemTopology::unmask_forces(const int atom_idx, const set<int> &substruct)
{
        for (auto &data : bondStretchData[atom_idx])
        { // get all forces involving this atom's idx
                if (substruct.count(data.idx1) && substruct.count(data.idx2))
                        bondStretch->setBondParameters(data.force_idx, data.idx1, data.idx2, data.length, data.k);
        }

        for (auto &data : bondBendData[atom_idx])
        { // get all forces involving this atom's idx
                if (substruct.count(data.idx1) && substruct.count(data.idx2) && substruct.count(data.idx3))
                        bondBend->setAngleParameters(data.force_idx, data.idx1, data.idx2, data.idx3, data.angle, data.k);
        }

        for (auto &data : bondTorsionData[atom_idx])
        { // get all forces involving this atom's idx
                if (substruct.count(data.idx1) && substruct.count(data.idx2) && substruct.count(data.idx3) && substruct.count(data.idx4))
                {

                        int idx1, idx2, idx3, idx4;
                        int periodicity;
                        double phase, k;

                        bondTorsion->getTorsionParameters(data.force_idx, idx1, idx2, idx3, idx4, periodicity, phase, k);
                        dbgmsg("particles of force " << data.force_idx);
                        dbgmsg(idx1 << " " << data.idx1);
                        dbgmsg(idx2 << " " << data.idx2);
                        dbgmsg(idx3 << " " << data.idx3);
                        dbgmsg(idx4 << " " << data.idx4);
                        dbgmsg(periodicity << " " << data.periodicity);

                        if (idx1 != data.idx1)
                                throw Error("die : particles changed");

                        if (idx2 != data.idx2)
                                throw Error("die : particles changed");

                        if (idx3 != data.idx3)
                                throw Error("die : particles changed");

                        if (idx4 != data.idx4)
                                throw Error("die : particles changed");

                        bondTorsion->setTorsionParameters(data.force_idx, data.idx1, data.idx2, data.idx3, data.idx4, data.periodicity, data.phase, data.k);
                }
        }
}

// Choose an Integrator for advancing time, and a Context connecting the
// System with the Integrator for simulation. Let the Context choose the
// best available Platform. Initialize the configuration from the default
// positions we collected above. Initial velocities will be zero.

void SystemTopology::init_integrator(SystemTopology::integrator_type type,
                                     const double step_size_in_ps,
                                     const double temperature_in_kevin,
                                     const double friction_in_per_ps)
{

        if (integrator != nullptr)
                throw Error("Integrator already initialized");

        __integrator_used = type;

        switch (__integrator_used)
        {
        case verlet:
                __thermostat_idx = system->addForce(new OpenMM::AndersenThermostat(temperature_in_kevin, friction_in_per_ps));
                // Fall through
        case none:
                integrator = new OpenMM::VerletIntegrator(step_size_in_ps);
                break;
        case langevin:
                integrator = new OpenMM::LangevinIntegrator(temperature_in_kevin, friction_in_per_ps, step_size_in_ps);
                break;
        case brownian:
                integrator = new OpenMM::BrownianIntegrator(temperature_in_kevin, friction_in_per_ps, step_size_in_ps);
                break;
        }

        map<string, string> properties;
        properties["CudaPrecision"] = "double";
        //properties["DisablePmeStream"] = "true";

        //int numParameters = forcefield->getNumPerParticleParameters();
        //std::cerr << "num particles per parameter " << numParameters << std::endl;
        //vector<double> parameters;
        //forcefield->getParticleParameters(0, parameters);
        //std::cerr << "parameter size " << parameters.size() << std::endl;

        //std::cerr << "step size " << step_size_in_ps << std::endl;
        OpenMM::Platform &platform = OpenMM::Platform::getPlatformByName("CUDA");
        context = new OpenMM::Context(*system, *integrator, platform, properties);
        //std::cerr << "Using OpenMM platform " << context->getPlatform().getName() << std::endl;

        dbgmsg("REMARK  Using OpenMM platform " << context->getPlatform().getName());
}

void SystemTopology::init_particles(Topology &topology)
{

        int warn = 0;

        system = new OpenMM::System();

        dbgmsg("topology.atoms.size() = " << topology.atoms.size());

        // Specify the atoms and their properties:
        //  (1) System needs to know the masses.
        //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
        //  (3) Collect default positions for initializing the simulation later.

        for (auto &patom : topology.atoms)
        {
                const Molib::Atom &atom = *patom;
                const int type = topology.get_type(atom);

                try
                {

                        const ForceField::AtomType &atype = __ffield->get_atom_type(type);
                        system->addParticle(atype.mass);

                        masses.push_back(atype.mass);
                        masked.push_back(false); // at the start unmask everything
                        dbgmsg("add particle type = " << type << " crd = " << atom.crd() << " mass = "
                                                      << atype.mass << " charge = " << atype.charge << " sigma = " << atype.sigma
                                                      << " epsilon = " << atype.epsilon << " representing atom = "
                                                      << atom << " at index = " << topology.get_index(atom));
                }
                catch (ParameterError &e)
                {
                        log_error << e.what() << " (" << ++warn << ")" << endl;
                }
        }

        if (warn > 0)
        {
                throw Error("die : missing parameters detected");
        }
}

void SystemTopology::update_thermostat(const double temperature_in_kevin,
                                       const double collision_frequency)
{

        switch (__integrator_used)
        {
        case verlet:
                dynamic_cast<OpenMM::AndersenThermostat &>(system->getForce(__thermostat_idx)).setDefaultTemperature(temperature_in_kevin);
                dynamic_cast<OpenMM::AndersenThermostat &>(system->getForce(__thermostat_idx)).setDefaultCollisionFrequency(collision_frequency);
                break;
        case langevin:
                dynamic_cast<OpenMM::LangevinIntegrator *>(integrator)->setTemperature(temperature_in_kevin);
                dynamic_cast<OpenMM::LangevinIntegrator *>(integrator)->setFriction(collision_frequency);
                break;
        case brownian:
                dynamic_cast<OpenMM::BrownianIntegrator *>(integrator)->setTemperature(temperature_in_kevin);
                dynamic_cast<OpenMM::BrownianIntegrator *>(integrator)->setFriction(collision_frequency);
                break;
        case none:
        default:
                throw Error("No integrator is being used!");
                break;
        }
}

void SystemTopology::clear_knowledge_based_force()
{
        //__kbforce->clearBonds();
        //__kbforce->updateParametersInContext(*context);
        forcefield->updateParametersInContext(*context);
}

void SystemTopology::update_knowledge_based_force(Topology &topology, const vector<OpenMM::Vec3> &positions, const double dist_cutoff)
{
        try
        {
                // we have to set new coordinates for atoms
                /* AtomPoint::UPVec atompoints;

                for (size_t i = 0; i < positions.size(); ++i)
                {
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

                for (auto &pa1 : atompoints)
                {
                        Molib::Atom &atom1 = pa1->get_atom();
                        visited.insert(&atom1);

                        for (auto &pa2 : grid.get_neighbors(pa1->crd(), dist_cutoff))
                        {
                                Molib::Atom &atom2 = pa2->get_atom();

                                if (!visited.count(&atom2))
                                {
                                        if (!topology.bonded_exclusions.count({&atom1, &atom2}) && !topology.bonded_exclusions.count({&atom2, &atom1}))
                                        {
                                                nonbondlist.push_back({&atom1, &atom2});
                                        }
                                }
                        }
                }

                dbgmsg("NONBOND LIST (KBFORCES) (SIZE = " << nonbondlist.size() << ") : " << endl
                                                          << nonbondlist);

                dbgmsg("kbforce step = " << setprecision(9) << fixed << __ffield->step);

                int warn = 0;

                //__kbforce->clearBonds();
                // init OpenMM's kbforce object
                for (auto &bond : nonbondlist)
                {
                        dbgmsg("next iteration in nonbondlist");
                        const Molib::Atom &atom1 = *bond.first;
                        const Molib::Atom &atom2 = *bond.second;
                        const int idx1 = topology.get_index(atom1);
                        const int idx2 = topology.get_index(atom2);

                        if (!masked[idx1] && !masked[idx2])
                        { // don't make the force if one or both atoms are masked
                                try
                                {
                                        //__kbforce->addBond (idx1, idx2, atom1.idatm_type(), atom2.idatm_type());
                                        dbgmsg("adding kbforce between idx1 = " << idx1 << " and idx2 = " << idx2);
                                }
                                catch (ParameterError &e)
                                {
                                        log_error << e.what() << " (" << ++warn << ")" << endl;
                                }
                        }
                } */

                dbgmsg("out of loop");
                //__kbforce->updateParametersInContext(*context);
                //forcefield->updateParametersInContext(*context);

                dbgmsg("after reinitialize");
                //context->setPositions(positions);
                dbgmsg("after setPositions");

                /* if (warn > 0)
                {
                        throw Error("die : missing parameters detected");
                }*/

                dbgmsg("exiting update_knowledge_based_force");
        }
        catch (std::exception &e)
        {
                log_error << e.what() << endl;
        }
}

void SystemTopology::init_physics_based_force(Topology &topology)
{

        int warn = 0;

        OpenMM::NonbondedForce *nonbond = new OpenMM::NonbondedForce();
        system->addForce(nonbond);

        for (auto &patom : topology.atoms)
        {
                const Molib::Atom &atom = *patom;
                const int type = topology.get_type(atom);

                try
                {
                        const ForceField::AtomType &atype = __ffield->get_atom_type(type);

                        nonbond->addParticle(atype.charge, atype.sigma, atype.epsilon);

                        dbgmsg("add particle type = " << type << " crd = " << atom.crd() << " mass = "
                                                      << atype.mass << " charge = " << atype.charge << " sigma = " << atype.sigma
                                                      << " epsilon = " << atype.epsilon << " representing atom = "
                                                      << atom << " at index = " << topology.get_index(atom));
                }
                catch (ParameterError &e)
                {
                        log_error << e.what() << " (" << ++warn << ")" << endl;
                }
        }

        vector<pair<int, int>> bondPairs;

        for (auto &bond : topology.bonds)
        {
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

        if (warn > 0)
        {
                throw Error("die : missing parameters detected");
        }
}

void SystemTopology::init_knowledge_based_force(Topology &topology)
{
        if (__kbforce_idx != -1)
                system->removeForce(__kbforce_idx);

        std::set<int> used_atom_types;

        forcefield = new OpenMM::CustomNonbondedForce("kbpot( idatm1, idatm2 , r / ffstep); ");
        forcefield->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        forcefield->setCutoffDistance(__ffield->kb_cutoff);
        //FIXME Fixup topology to contain all the atom types used in the system (thus no need to initialize per ligand/protein)
        std::map<int, int> __idatm_to_internal;
        std::map<int, int> __internal_to_idatm;
        int num_types = 0;
        for (const auto &atom : topology.atoms)
        {
                used_atom_types.insert(atom->idatm_type());
                if (!__idatm_to_internal.count(atom->idatm_type()))
                {
                        __idatm_to_internal[atom->idatm_type()] = num_types;
                        __internal_to_idatm[num_types] = atom->idatm_type();
                        num_types++;
                }

                forcefield->addParticle({__idatm_to_internal[atom->idatm_type()]});
        }
        forcefield->addGlobalParameter("ffstep", __ffield->step);
        forcefield->addPerParticleParameter("idatm");
        /* Disable KBForce */

        // __kbforce = new KBPlugin::KBForce(used_atom_types.size(), __ffield->step, __ffield->kb_cutoff);

        vector<double> table;
        int xsize = 0, ysize = 0;

        for (int i = 0; i < __idatm_to_internal.size(); i++)
        {
                for (int j = 0; j < __idatm_to_internal.size(); j++)
                {
                        try
                        {
                                const ForceField::KBType kb = __ffield->get_kb_force_type(__internal_to_idatm[i], __internal_to_idatm[j]);

                                //std::cerr << "size " << kb.potential.size() << std::endl;
                                if (ysize == 0)
                                        ysize = kb.potential.size();

                                if (ysize != kb.potential.size())
                                        throw Error("mismatching potential size");

                                xsize++;

                                for (size_t i = 0; i < kb.potential.size(); i++)
                                        table.push_back(kb.potential[i]);
                                //__kbforce->addBondType( type1, type2, kb.potential, kb.derivative);
                        }
                        catch (ParameterError &e)
                        {
                                cerr << e.what() << endl;
                                cerr << "This is normal for atom types only present in the ligand" << endl;

                                //size_t number_of_steps = static_cast<size_t>(__ffield->kb_cutoff / __ffield->step);
                                // __kbforce->addBondType( type1, type2, std::vector<double>(number_of_steps, 0.0),
                                //    std::vector<double>(number_of_steps, 0.0) );addForce
                        }
                }
        }
        try
        {
                //std::cerr << "x y table " << xsize << " " << ysize << " " << table.size() << std::endl;
                //forcefield->addGlobalParameter("ntypes", __idatm_to_internal.size());

                forcefield->addTabulatedFunction("kbpot", new OpenMM::Discrete3DFunction(__idatm_to_internal.size(), __idatm_to_internal.size(), ysize, table));

                vector<pair<int, int>> bondPairs;

                for (auto &bond : topology.bonds)
                {
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

                forcefield->createExclusionsFromBonds(bondPairs, 4);
                auto num_exclusions = forcefield->getNumExclusions();
                //std::cerr << "num exlcusions " << num_exclusions << std::endl;
                /*for (auto i = 0; i < forcefield->getNumExclusions(); i++)
                {
                        int p1 = -1, p2 = -1;
                        forcefield->getExclusionParticles(i, p1, p2);
                        std::cerr << "exclusions p1 " << p1 << " p2 " << p2 << std::endl;
                }*/
        }
        catch (ParameterError &e)
        {
                cerr << e.what() << endl;
                cerr << "Exiting" << endl;
                exit(0);
        }

        __kbforce_idx = system->addForce(forcefield);
}

void SystemTopology::retype_amber_protein_atom_to_gaff(const Molib::Atom &atom, int &type)
{

        // Only retype protein atoms to gaff atoms
        if (atom.br().rest() != Molib::Residue::protein)
        {
                // Its the N in the peptide bond needs to be N, not n3....
                if (atom.element() == Molib::Element::N)
                {
                        type = __ffield->gaff_name_to_type.at("n");
                }
        }

        if (atom.element() == Molib::Element::C)
        {
                // C2 == 18 aka the carbon is in the carbonyl
                if (atom.idatm_type() == 18)
                        type = __ffield->gaff_name_to_type.at("c");
                // Otherwise it's the alpha carbon
                else
                        type = __ffield->gaff_name_to_type.at("c3");
        }

        // Its the N in the peptide bond
        if (atom.element() == Molib::Element::N)
        {
                type = __ffield->gaff_name_to_type.at("n");
        }

        // Its the O in the peptide bond
        if (atom.element() == Molib::Element::O)
        {
                type = __ffield->gaff_name_to_type.at("o");
        }
}

void SystemTopology::init_bonded(Topology &topology, const bool use_constraints)
{

        int warn = 0;

        bondStretchData.resize(topology.atoms.size());
        bondBendData.resize(topology.atoms.size());
        bondTorsionData.resize(topology.atoms.size());

        bondStretch = new OpenMM::HarmonicBondForce();
        bondBend = new OpenMM::HarmonicAngleForce();
        bondTorsion = new OpenMM::PeriodicTorsionForce();

        system->addForce(bondStretch);
        system->addForce(bondBend);
        system->addForce(bondTorsion);

        dbgmsg("initializing openmm");

        // Process the bonds:
        //  (1) If we're using constraints, tell System about constrainable bonds;
        //      otherwise, tell HarmonicBondForce the bond stretch parameters
        //      (tricky units!).
        //  (2) Create a list of bonds for generating nonbond exclusions.

        int force_idx = 0;

        for (auto &bond : topology.bonds)
        {
                dbgmsg("checkpoint0");

                const Molib::Atom &atom1 = *bond.first;
                dbgmsg(atom1);
                const Molib::Atom &atom2 = *bond.second;
                dbgmsg(atom2);
                const int idx1 = topology.get_index(atom1);
                dbgmsg(idx1);
                const int idx2 = topology.get_index(atom2);
                dbgmsg(idx2);
                int type1 = topology.get_type(atom1);
                dbgmsg(type1);
                int type2 = topology.get_type(atom2);
                dbgmsg(type2);

                // Check for modified residues
                int number_protein = (atom1.br().rest() == Molib::Residue::protein) + (atom2.br().rest() == Molib::Residue::protein);

                int number_not_set = (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???");

                // It should both protein or non-protein, otherwise we fix it
                if (number_protein == 1 && number_not_set == 1)
                {
                        retype_amber_protein_atom_to_gaff(atom1, type1);
                        retype_amber_protein_atom_to_gaff(atom2, type2);
                }

                ForceField::BondType btype;

                try
                {
                        btype = __ffield->get_bond_type(type1, type2);
                }
                catch (ParameterError &e)
                {
                        log_warning << e.what() << " (WARNINGS ARE NOT INCREASED) (using default parameters for this bond)" << endl;
                        // if everything else fails just constrain at something reasonable
                        btype = ForceField::BondType{atom1.get_bond(atom2).length(), 250000, false};
                }

                if (use_constraints && btype.can_constrain)
                { // Should we constrain C-H bonds?
                        system->addConstraint(idx1, idx2, btype.length);
                }
                else
                {
                        // Note factor of 2 for stiffness below because Amber specifies the constant
                        // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
                        // it as used in the force term kx, with energy kx^2/2.
                        bondStretch->addBond(idx1, idx2, btype.length, btype.k);
                        dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1 << " idx2 = " << idx2
                                              << " bond length = " << btype.length << " k = " << btype.k);
                        bondStretchData[idx1].push_back(ForceData{force_idx, idx1, idx2, 0, 0, btype.length, 0, 0, 0, btype.k});
                        bondStretchData[idx2].push_back(ForceData{force_idx, idx1, idx2, 0, 0, btype.length, 0, 0, 0, btype.k});
                        ++force_idx;
                }
        }

        dbgmsg("checkpoint3");

        force_idx = 0;

        // Create the 1-2-3 bond angle harmonic terms.
        for (auto &angle : topology.angles)
        {
                dbgmsg("checkpoint4");
                const Molib::Atom &atom1 = *get<0>(angle);
                const Molib::Atom &atom2 = *get<1>(angle);
                const Molib::Atom &atom3 = *get<2>(angle);
                const int idx1 = topology.get_index(atom1);
                const int idx2 = topology.get_index(atom2);
                dbgmsg("checkpoint5");
                const int idx3 = topology.get_index(atom3);
                int type1 = topology.get_type(atom1);
                int type2 = topology.get_type(atom2);
                int type3 = topology.get_type(atom3);
                dbgmsg("checkpoint6");

                // Check for modified residues
                int number_protein = (atom1.br().rest() == Molib::Residue::protein) + (atom2.br().rest() == Molib::Residue::protein) + (atom3.br().rest() == Molib::Residue::protein);

                int number_not_set = (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???") +
                                     (atom3.gaff_type() != "???");

                // It should all be protein or non-protein, otherwise we fix it
                if ((number_protein == 1 || number_protein == 2) && (number_not_set == 1 || number_not_set == 2))
                {
                        retype_amber_protein_atom_to_gaff(atom1, type1);
                        retype_amber_protein_atom_to_gaff(atom2, type2);
                        retype_amber_protein_atom_to_gaff(atom3, type3);
                }

                ForceField::AngleType atype;

                try
                {
                        dbgmsg("determining angle type between atoms : " << endl
                                                                         << atom1 << endl
                                                                         << atom2 << endl
                                                                         << atom3);
                        atype = __ffield->get_angle_type(type1, type2, type3);
                }
                catch (ParameterError &e)
                {
                        log_warning << e.what() << " (WARNINGS ARE NOT INCREASED) (using default parameters for this angle)" << endl;
                        // if everything else fails just constrain at something reasonable
                        atype = ForceField::AngleType{Geom3D::angle(atom1.crd(), atom2.crd(), atom3.crd()), 500};
                }

                dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1 << " idx2 = " << idx2 << " idx3 = " << idx3
                                      << " angle = " << atype.angle << " k = " << atype.k);
                bondBend->addAngle(idx1, idx2, idx3, atype.angle, atype.k);
                bondBendData[idx1].push_back(ForceData{force_idx, idx1, idx2, idx3, 0, 0, atype.angle, 0, 0, atype.k});
                bondBendData[idx2].push_back(ForceData{force_idx, idx1, idx2, idx3, 0, 0, atype.angle, 0, 0, atype.k});
                bondBendData[idx3].push_back(ForceData{force_idx, idx1, idx2, idx3, 0, 0, atype.angle, 0, 0, atype.k});

                ++force_idx;

                dbgmsg("checkpoint8");
        }

        force_idx = 0;

        // Create the 1-2-3-4 bond torsion (dihedral) terms.
        for (auto &dihedral : topology.dihedrals)
        {
                dbgmsg("checkpoint9");
                const Molib::Atom &atom1 = *get<0>(dihedral);
                const Molib::Atom &atom2 = *get<1>(dihedral);
                const Molib::Atom &atom3 = *get<2>(dihedral);
                const Molib::Atom &atom4 = *get<3>(dihedral);
                const int idx1 = topology.get_index(atom1);
                const int idx2 = topology.get_index(atom2);
                const int idx3 = topology.get_index(atom3);
                const int idx4 = topology.get_index(atom4);
                int type1 = topology.get_type(atom1);
                int type2 = topology.get_type(atom2);
                int type3 = topology.get_type(atom3);
                int type4 = topology.get_type(atom4);

                // Check for modified residues
                int number_protein = (atom1.br().rest() == Molib::Residue::protein) + (atom2.br().rest() == Molib::Residue::protein) + (atom3.br().rest() == Molib::Residue::protein) + (atom4.br().rest() == Molib::Residue::protein);

                int number_not_set = (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???") +
                                     (atom3.gaff_type() != "???") + (atom4.gaff_type() != "???");

                // It should all be protein or non-protein, otherwise we fix it
                if ((number_protein == 1 || number_protein == 2 || number_protein == 3) && (number_not_set == 1 || number_not_set == 2 || number_not_set == 3))
                {
                        retype_amber_protein_atom_to_gaff(atom1, type1);
                        retype_amber_protein_atom_to_gaff(atom2, type2);
                        retype_amber_protein_atom_to_gaff(atom3, type3);
                        retype_amber_protein_atom_to_gaff(atom4, type4);
                }

                try
                {
                        const ForceField::TorsionTypeVec &v_ttype = __ffield->get_dihedral_type(type1, type2, type3, type4); // cannot make it const ??
                        for (auto &ttype : v_ttype)
                        {
                                dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1 << " idx2 = " << idx2
                                                      << " idx3 = " << idx3 << " idx4 = " << idx4 << " periodicity = " << ttype.periodicity
                                                      << " phase = " << ttype.phase << " k = " << ttype.k);
                                bondTorsion->addTorsion(idx1, idx2, idx3, idx4, ttype.periodicity, ttype.phase, ttype.k);
                                bondTorsionData[idx1].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});
                                bondTorsionData[idx2].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});
                                bondTorsionData[idx3].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});
                                bondTorsionData[idx4].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});

                                ++force_idx;
                        }
                }
                catch (ParameterError &e)
                {
                        log_warning << e.what() << " (WARNINGS ARE NOT INCREASED) (using default parameters for this dihedral)" << endl;
                        // if everything else fails just constrain at something reasonable
                        //cout <<  Geom3D::dihedral(atom1.crd(), atom2.crd(), atom3.crd(), atom4.crd()) << endl;
                        //ForceField::TorsionType vtype = ForceField::TorsionType {1, Geom3D::dihedral(atom1.crd(), atom2.crd(), atom3.crd(), atom4.crd()), 0};
                        //bondTorsion->addTorsion (idx1, idx2, idx3, idx4, vtype.periodicity, vtype.phase, vtype.k);
                        //bondTorsionData[idx1].push_back (ForceData {force_idx, idx1, idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase, vtype.k});
                        //bondTorsionData[idx2].push_back (ForceData {force_idx, idx1, idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase, vtype.k});
                        //bondTorsionData[idx3].push_back (ForceData {force_idx, idx1, idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase, vtype.k});
                        //bondTorsionData[idx4].push_back (ForceData {force_idx, idx1, idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase, vtype.k});

                        //++force_idx;
                }
        }

        // Create the 1-2-3-4 improper terms where 3 is the central atom
        for (auto &dihedral : topology.impropers)
        {
                dbgmsg("checkpoint10");
                const Molib::Atom &atom1 = *get<0>(dihedral);
                const Molib::Atom &atom2 = *get<1>(dihedral);
                const Molib::Atom &atom3 = *get<2>(dihedral);
                const Molib::Atom &atom4 = *get<3>(dihedral);
                const int idx1 = topology.get_index(atom1);
                const int idx2 = topology.get_index(atom2);
                const int idx3 = topology.get_index(atom3);
                const int idx4 = topology.get_index(atom4);
                int type1 = topology.get_type(atom1);
                int type2 = topology.get_type(atom2);
                int type3 = topology.get_type(atom3);
                int type4 = topology.get_type(atom4);

                // Check for modified residues
                int number_protein = (atom1.br().rest() == Molib::Residue::protein) + (atom2.br().rest() == Molib::Residue::protein) + (atom3.br().rest() == Molib::Residue::protein) + (atom4.br().rest() == Molib::Residue::protein);

                int number_not_set = (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???") +
                                     (atom3.gaff_type() != "???") + (atom4.gaff_type() != "???");

                // It should all be protein or non-protein, otherwise we fix it
                if ((number_protein == 1 || number_protein == 2 || number_protein == 3) && (number_not_set == 1 || number_not_set == 2 || number_not_set == 3))
                {
                        retype_amber_protein_atom_to_gaff(atom1, type1);
                        retype_amber_protein_atom_to_gaff(atom2, type2);
                        retype_amber_protein_atom_to_gaff(atom3, type3);
                        retype_amber_protein_atom_to_gaff(atom4, type4);
                }

                try
                {
                        const ForceField::TorsionTypeVec &v_ttype =
                            __ffield->get_improper_type(type1, type2, type3, type4); // cannot make it const ??

                        for (auto &ttype : v_ttype)
                        {
                                dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1 << " idx2 = " << idx2
                                                      << " idx3 = " << idx3 << " idx4 = " << idx4);
                                bondTorsion->addTorsion(idx1, idx2, idx3, idx4, ttype.periodicity, ttype.phase, ttype.k);
                                bondTorsionData[idx1].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});
                                bondTorsionData[idx2].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});
                                bondTorsionData[idx3].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});
                                bondTorsionData[idx4].push_back(ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0, ttype.periodicity, ttype.phase, ttype.k});

                                ++force_idx;
                        }
                }
                catch (ParameterError &e)
                {
                        dbgmsg(e.what() << " (WARNINGS ARE NOT INCREASED)");
                }
        }

        if (warn > 0)
        {
                throw Error("die : missing parameters detected");
        }
}

void SystemTopology::init_positions(const Geom3D::Point::Vec &crds)
{
        dbgmsg("entering init_positions crds.size() = " << crds.size());
        vector<OpenMM::Vec3> positions_in_nm;

        for (auto &crd : crds)
        {
                positions_in_nm.push_back(OpenMM::Vec3(crd.x() * OpenMM::NmPerAngstrom,
                                                       crd.y() * OpenMM::NmPerAngstrom, crd.z() * OpenMM::NmPerAngstrom));
        }

        context->setPositions(positions_in_nm);

        dbgmsg("exiting init_positions");
}

vector<OpenMM::Vec3> SystemTopology::get_positions_in_nm()
{
        return context->getState(OpenMM::State::Positions).getPositions();
}

vector<OpenMM::Vec3> SystemTopology::get_forces()
{
        return context->getState(OpenMM::State::Forces).getForces();
}

double SystemTopology::get_potential_energy()
{
        return context->getState(OpenMM::State::Energy).getPotentialEnergy();
}

void SystemTopology::minimize(const double tolerance, const double max_iterations)
{
        OpenMM::LocalEnergyMinimizer::minimize(*context, tolerance, max_iterations);
}

void SystemTopology::dynamics(const int steps)
{
        if (__integrator_used == integrator_type::verlet && __thermostat_idx == -1)
                log_warning << "No thermostat set, performing NVE dynamics" << endl;

        integrator->step(steps);
}
};
