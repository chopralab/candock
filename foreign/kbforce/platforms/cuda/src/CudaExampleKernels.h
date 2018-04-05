#ifndef CUDA_EXAMPLE_KERNELS_H_
#define CUDA_EXAMPLE_KERNELS_H_

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
#include "../../../openmmapi/include/KBKernels.h"
#include "KBKernels.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaArray.h"
#include <vector>


namespace KBPlugin {

/**
 * This kernel is invoked by KBForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcExampleForceKernel : public CalcKBForceKernel {
public:
    CudaCalcExampleForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::CudaContext& cu, const OpenMM::System& system) :
            CalcKBForceKernel(name, platform), hasInitializedKernel(false), cu(cu), system(system) {
    }
    ~CudaCalcExampleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the KBForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const KBForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ExampleForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const KBForce& force);
private:

    void copyBonds( const KBForce& force);    

    bool hasInitializedKernel;
    OpenMM::CudaContext& cu;
    const OpenMM::System& system;
    //OpenMM::CudaArray* params;
    OpenMM::CudaNonbondedUtilities* nonbonded;
    
    OpenMM::CudaArray *dev_potential_lookup, *dev_derivative_lookup,
        *dev_particle1, *dev_particle2, *dev_type1, *dev_type2;

    OpenMM::CudaArray* tempPosq;
    OpenMM::CudaArray* tempForces;
        
   // int numBonds;
    const double*  potential_lookup;
    const double* derivative_lookup;
    std::vector<int> particle1, particle2;
    std::vector<int> type1, type2;
    int    number_of_types;
    int    number_of_interactions;
    double step;
    double cutoff;
    size_t number_of_steps;
    int numBonds;


};

} // namespace KBPlugin

#endif /*CUDA_EXAMPLE_KERNELS_H_*/
