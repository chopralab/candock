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

#include "ReferenceKBKernels.h"
#include "KBForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include "helper/debug.hpp"
#include <iostream>

using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

namespace KBPlugin {

    void ReferenceCalcKBForceKernel::initialize(const System& system, const KBForce& force) {
        // Initialize bond parameters.
        dbgmsg("initializing bond parameters");

        int numBonds = force.getNumBonds();
        particle1.resize(numBonds);
        particle2.resize(numBonds);
        //~ length.resize(numBonds);
        //~ k.resize(numBonds);
        potential.resize(numBonds);
        derivative.resize(numBonds);
        for (int i = 0; i < numBonds; i++)
        //~ force.getBondParameters(i, particle1[i], particle2[i], length[i], k[i]);
        force.getBondParameters(i, particle1[i], particle2[i], potential[i], derivative[i]);
        // initialize global parameters.
        step = force.getStep();
    }

    double ReferenceCalcKBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
        dbgmsg("starting energy calculation (in execute) ");
        vector<RealVec>& pos = extractPositions(context);
        vector<RealVec>& force = extractForces(context);
        int numBonds = particle1.size();
        double energy = 0;

        dbgmsg("numBonds = " << numBonds);
        // Compute the interactions.
    
        for (int i = 0; i < numBonds; i++) {
            int p1 = particle1[i];
            int p2 = particle2[i];
            RealVec delta = pos[p1] - pos[p2];
            RealOpenMM r2 = delta.dot(delta);
            RealOpenMM r = sqrt(r2);
            size_t dist = static_cast<size_t> (floor(r / step));
            if (dist >= (*potential[i]).size()) continue; // effectively add zero to energy and force
            energy += (*potential[i])[dist];

            RealOpenMM dEdR = (*derivative[i])[dist];
            dEdR = (r > 0) ? (dEdR/r) : 0;
            force[p1] -= delta*dEdR;
            force[p2] += delta*dEdR;

            dbgmsg("kb particle1 = " << p1 << " particle2 = " << p2 
                << " distance = " << r << " index for potential = " << dist 
                << " potential = " << (*potential[i])[dist] 
                << " derivative = " << (*derivative[i])[dist] 
                << " size of potential = " << (*potential[i]).size()
                << " force[" << p1 << "] = " << force[p1]
                << " force[" << p2 << "] = " << force[p2]);

        }
        dbgmsg("knowledge-based energy is calculated as : energy = " << setprecision(20) << energy);
        return energy;
    }

    void ReferenceCalcKBForceKernel::copyParametersToContext(ContextImpl& context, const KBForce& force) {
        if ( static_cast<size_t>(force.getNumBonds()) != particle1.size())
            throw OpenMMException("updateParametersInContext: The number of KB bonds has changed");
        for (int i = 0; i < force.getNumBonds(); i++) {
            int p1, p2;
            //~ force.getBondParameters(i, p1, p2, length[i], k[i]);
            force.getBondParameters(i, p1, p2, potential[i], derivative[i]);
            if (p1 != particle1[i] || p2 != particle2[i])
                throw OpenMMException("updateParametersInContext: A particle index has changed");
        }
    }
}