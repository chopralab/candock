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
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "openmm/reference/RealVec.h"


#include <time.h>

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
    
    if(dev_particle1)
        delete dev_particle1;
    if(dev_particle2)
        delete dev_particle2;
    if(dev_type1)
        delete dev_type1;
    if(dev_type2)
        delete dev_type2;
}

void CudaCalcExampleForceKernel::initialize(const System &system, const KBForce &force)
{
    
    cu.setAsCurrent();
    
    potential_lookup  = force.getPotentialLookup();
    derivative_lookup = force.getDerivativeLookup();
    numBonds = force.getNumBonds();
    number_of_types = force.getNumTypes();
    number_of_interactions = force.getNumInteractions();
    step = force.getStep();
    cutoff = force.getCutoff();
    number_of_steps = force.getNumSteps();
    
    //cout << "Number of steps " << number_of_steps << " Number of Interactions " << number_of_interactions << endl;
    dev_potential_lookup = CudaArray::create<double>(cu, number_of_steps * number_of_interactions, "dev_potential_lookup");
    dev_derivative_lookup = CudaArray::create<double>(cu, number_of_steps * number_of_interactions, "dev_derivative_lookup");
    cudaDeviceSynchronize();
    dev_potential_lookup->upload(potential_lookup);
    dev_derivative_lookup->upload(derivative_lookup);
    
    
    cu.addForce(new CudaExampleForceInfo(force));
    
    
}



double CudaCalcExampleForceKernel::execute(ContextImpl &context, bool includeForces, bool includeEnergy) 
{ 
    
    //Handle the initialize step, which calls execute, but has 0 bonds.
    if (numBonds <= 0) {
        //cout << "Energy 0" << endl;
        return 0.0;
    }
    
    
    double *energy = new double[1];
    energy[0] = 0.0;
    OpenMM::CudaArray *dev_energy = CudaArray::create<double>(cu, 1, "dev_energy");
    dev_energy->upload(energy);
    
    //Defines enables you to access vars in the kernel
    map<string, string> defines;
    defines["num_bonds"] = cu.intToString(numBonds);
    defines["step"] = cu.doubleToString(step);
    defines["number_of_types"] = cu.intToString(number_of_types);
    defines["number_of_interactions"] = cu.intToString(number_of_interactions);
    defines["cutoff"] = cu.doubleToString(cutoff);
    defines["number_of_steps"] = cu.intToString(static_cast<int>(number_of_steps));
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["step"] = cu.doubleToString(step);
    
    cu.getUseDoublePrecision();
    CUmodule module = cu.createModule(CudaExampleKernelSources::exampleForce, defines, "");
    CUfunction kernel = cu.getKernel(module, "calcNonBonded");

    void* prepareArgs[] = {
        &cu.getPosq().getDevicePointer(),
        &cu.getForce().getDevicePointer(),
        &dev_particle1->getDevicePointer(), 
        &dev_particle2->getDevicePointer(), 
        &dev_type1->getDevicePointer(), 
        &dev_type2->getDevicePointer(), 
        &dev_potential_lookup->getDevicePointer(), 
        &dev_derivative_lookup->getDevicePointer(),
        &dev_energy->getDevicePointer(),
    };
    cudaDeviceSynchronize();
    
    cu.executeKernel(kernel, prepareArgs, static_cast<unsigned int>(numBonds));
    
    cudaDeviceSynchronize();
    
    dev_energy->download(energy);
    
    cudaDeviceSynchronize();
   // cout << "Energy " << energy[0] << endl;
    return energy[0];
}

void CudaCalcExampleForceKernel::copyParametersToContext(ContextImpl &context, const KBForce &force)
{
    cu.setAsCurrent();
    particle1.clear();
    particle2.clear();
    type1.clear();
    type2.clear();
    
    numBonds = force.getNumBonds();
    if (numBonds <= 0) {
        return;
    }
    
    dev_particle1 = CudaArray::create<int>(cu, (numBonds), "dev_particle1");
    dev_particle2 = CudaArray::create<int>(cu, (numBonds), "dev_particle2");
    dev_type1 = CudaArray::create<int>(cu, (numBonds), "dev_type1");
    dev_type2 = CudaArray::create<int>(cu, (numBonds), "dev_type2");
    
    
    particle1.resize(numBonds);
    particle2.resize(numBonds);
    type1.resize(numBonds);
    type2.resize(numBonds);
    
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, particle1[i], particle2[i], type1[i], type2[i]);
    
    dev_particle2->upload(particle2);
    dev_particle1->upload(particle1);
    
    dev_type1->upload(type1);
    dev_type2->upload(type2);
    cudaDeviceSynchronize();
    
    cu.invalidateMolecules();
    
}

