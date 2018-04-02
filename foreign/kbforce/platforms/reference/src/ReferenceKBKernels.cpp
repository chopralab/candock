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
#include "openmm/reference/ReferencePlatform.h"
#include <iostream>
#include <algorithm>

using namespace KBPlugin;
using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

void ReferenceCalcKBForceKernel::initialize(const System& system, const KBForce& force) {
    copyBonds(force);
}

double ReferenceCalcKBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
 
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& force = extractForces(context);
    size_t numBonds = particle1.size();
    double energy = 0;

    // Compute the interactions.
    
    for (size_t i = 0; i < numBonds; i++) {
        int p1 = particle1[i];
        int p2 = particle2[i];
        Vec3 delta = pos[p1] - pos[p2];
        double r2 = delta.dot(delta);
        double r = sqrt(r2);
        int dist = std::floor(r / step);

        const int internal_min = std::min(type1[i], type2[i]);
        const int internal_max = std::max(type1[i], type2[i]);
        const int offset = (internal_max + internal_min*(number_of_types-1)-internal_min*(internal_min-1)/2) * number_of_steps;

        if (dist >= number_of_steps) continue; // effectively add zero to energy and force
        energy += potential_lookup[offset + dist];

        double dEdR = derivative_lookup[offset + dist];
        dEdR = (r > 0) ? (dEdR/r) : 0;
        force[p1] -= delta*dEdR;
        force[p2] += delta*dEdR;

    }
    return energy;
}

void ReferenceCalcKBForceKernel::copyParametersToContext(ContextImpl& context, const KBForce& force) {

    particle1.clear();
    particle2.clear();
    type1.clear();
    type2.clear();

    copyBonds(force);
}

void KBPlugin::ReferenceCalcKBForceKernel::copyBonds(const KBPlugin::KBForce& force) {

    potential_lookup = force.getPotentialLookup();
    derivative_lookup= force.getDerivativeLookup();
        
    int numBonds = force.getNumBonds();

    particle1.resize(numBonds);
    particle2.resize(numBonds);
    type1.resize(numBonds);
    type2.resize(numBonds);

    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, particle1[i], particle2[i], type1[i], type2[i]);

    number_of_types = force.getNumTypes();
    number_of_interactions = force.getNumInteractions();
    step = force.getStep();
    cutoff = force.getCutoff();
    number_of_steps = force.getNumSteps();
}
