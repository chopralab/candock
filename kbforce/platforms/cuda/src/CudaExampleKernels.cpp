#include <OpenMM.h>
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */
#include <iostream>
#include "CudaExampleKernels.h"
//#include <../openmm/include/KBForce.h>
#include <../cuda-8.0/include/vector_functions.hpp>
#include "CudaExampleKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaNonbondedUtilities.h"
#include "openmm/reference/RealVec.h"

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <numeric>
#include <time.h>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <ratio>
#include <math.h>
#include <unistd.h>

using namespace KBPlugin;
using namespace OpenMM;
using namespace std;

class CudaExampleForceInfo : public CudaForceInfo
{
  public:
    CudaExampleForceInfo(const KBForce &force) : force(force)
    {
    }

  private:
    const KBForce &force;
};

CudaCalcExampleForceKernel::~CudaCalcExampleForceKernel()
{
    cu.setAsCurrent();

}

/*
 **************************************************************************************
 *
 *  initialize - This is the first step run inside this OpenMM plugin.
 * 
 *  This function is designed to set up all the necessary variables and create the
 *  lookup tables. These lookup tables are used in the execute function for determining
 *  the force and energy.
 * 
 **************************************************************************************
 */


void CudaCalcExampleForceKernel::initialize(const System &system, const KBForce &force)
{
    cerr << "Entered Initialize" << endl;
    
    cu.setAsCurrent();

    /* Retrieve Information From OpenMM  */
    potential_lookup = force.getPotentialLookup();
    derivative_lookup = force.getDerivativeLookup();
    number_of_interactions = force.getNumInteractions();
    number_of_steps = force.getNumSteps();
    forceGroup = force.getForceGroup();
    this->run_non_bonded_interactions = false;
    cerr << "num atoms" << cu.getNumAtoms() << endl;
    cerr << forceGroup << endl;

    /* Create potential look up table and derivative lookup table */
    dev_potential_lookup = CudaArray::create<double>(cu, number_of_steps * number_of_interactions, "dev_potential_lookup");
    dev_derivative_lookup = CudaArray::create<double>(cu, number_of_steps * number_of_interactions, "dev_derivative_lookup");

    dev_potential_lookup->upload(potential_lookup);
    dev_derivative_lookup->upload(derivative_lookup);
    
    nonbonded = &cu.getNonbondedUtilities();
   
    nonbonded->addArgument(CudaNonbondedUtilities::ParameterInfo("dev_potential_lookup", "double", number_of_steps * number_of_interactions,
        sizeof(double), dev_potential_lookup->getDevicePointer()));
    nonbonded->addArgument(CudaNonbondedUtilities::ParameterInfo("dev_derivative_lookup", "double", number_of_steps * number_of_interactions,
        sizeof(double), dev_derivative_lookup->getDevicePointer()));
}

/*
 **************************************************************************************
 *
 *  execute - calculate nonbonded forces and reduce the energy.
 * 
 **************************************************************************************
 */


double CudaCalcExampleForceKernel::execute(ContextImpl &context, bool includeForces, bool includeEnergy)
{

    cerr << "num atoms" << cu.getNumAtoms() << endl;
    cerr << "Entered Execute" << endl;
 
    if(!run_non_bonded_interactions) return 0.0;


    cu.setAsCurrent();

    if (!hasInitializedNonbonded) {
        hasInitializedNonbonded = true;
        nonbonded->initialize(system);
    }

    /* Calcualte nonbonded forces */
    nonbonded->computeInteractions(forceGroup, 1, 1);

    cerr << "calculated nonbonded forces" << endl;
    double energy = cu.reduceEnergy();
    cerr << "energy " << energy << endl;


    return energy;
}


/*
 **************************************************************************************
 *
 *   copyParametersToContext - Update the neighbor lists
 * 
 **************************************************************************************
 */

void CudaCalcExampleForceKernel::copyParametersToContext(ContextImpl &context, const KBForce &force, const double dist_cutoff)
{

    cerr << "num atoms" << cu.getNumAtoms() << endl;
    cerr << "Entered Copy Parameters" << endl;

    if(!run_non_bonded_interactions) return;


    cu.setAsCurrent();

     if (!hasInitializedNonbonded) {
        hasInitializedNonbonded = true;
        nonbonded->initialize(system);
    }

    cerr << "Start calculating neighbor lists" << endl;
    /* Compute Neighbor Lists */
    nonbonded->prepareInteractions(force.getForceGroup());

    cerr << "Finished calculating neighbor lists" << endl;

    cu.invalidateMolecules();

    cerr << "leave Copy Parameters" << endl;
}


/*
 **************************************************************************************
 *
 *  copyTypes - 
 * 
 *  Runs once a few steps after initialize is called, but before execute or 
 *  copyParametersToContext is called. This function adds the bonded_exclusions_matrix
 *  to KBForce. It also creates the interaction.
 * 
 **************************************************************************************
 */


void CudaCalcExampleForceKernel::copyTypes(ContextImpl &context, int *atoms, int num_atoms, vector<vector<int>> bonded_exclusions_matrix, const KBForce& force)
{
    cerr << "Entered Copy Types" << endl;
    
    cerr << "num atoms" << cu.getNumAtoms() << endl;


    cu.setAsCurrent();

    /* Convert atoms to internal version used in OpenMM */
    /*vector<int> atoms_internal;
    for (size_t i = 0; i < num_atoms; i++)
        atoms_internal.push_back(force.convert(atoms[i]));
    
    dev_atoms = CudaArray::create<int>(cu, num_atoms, "dev_atoms");
    dev_atoms->upload(atoms_internal);
    

    nonbonded->addArgument(CudaNonbondedUtilities::ParameterInfo("dev_atoms", "int", num_atoms,
        sizeof(double), dev_atoms->getDevicePointer()));
    */

    map<string, string> defines;

    defines["step"] = cu.doubleToString(step);
    defines["number_of_types"] = cu.intToString(number_of_types);
    defines["number_of_interactions"] = cu.intToString(number_of_interactions);
    defines["number_of_steps"] = cu.intToString(static_cast<int>(number_of_steps));

    nonbonded->addInteraction(true, false, true, 0.06, bonded_exclusions_matrix, cu.replaceStrings(CudaExampleKernelSources::exampleForce, defines), force.getForceGroup());

    cu.addForce(new CudaExampleForceInfo(force));
}

void CudaCalcExampleForceKernel::calculate_non_bonded_forces()
{
    this->run_non_bonded_interactions = true;
}