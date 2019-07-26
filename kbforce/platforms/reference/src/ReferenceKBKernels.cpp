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
#include <iostream>
#include <algorithm>

#include <numeric>
#include <time.h>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <ratio>

using namespace KBPlugin;
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

void ReferenceCalcKBForceKernel::initialize(const System& system, const KBForce& force) {
    copyBonds(force);
}

double ReferenceCalcKBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //auto t5 = std::chrono::high_resolution_clock::now();   
    
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    size_t numBonds = particle1.size();

    double energy = 0;

    // Compute the interactions.
    for (size_t i = 0; i < numBonds; i++) {
        int p1 = particle1[i];
        int p2 = particle2[i];
        RealVec delta = pos[p1] - pos[p2];
        RealOpenMM r2 = delta.dot(delta);
        RealOpenMM r = sqrt(r2);
        int dist = std::floor(r / step);

        const int internal_min = std::min(type1[i], type2[i]);
        const int internal_max = std::max(type1[i], type2[i]);
        const int offset = (internal_max + internal_min*(number_of_types-1)-internal_min*(internal_min-1)/2) * number_of_steps;

        //fprintf(stdout, "i %d, pos2.x %lf, pos2.y %lf, pos2.z %lf, pos1.x %lf, pos1.y %lf, pos1.z %lf, offset %d, dist %d, p1 %d, p2 %d, energy %lf, r %lf, step %lf \n" , i, pos[p2][0], pos[p2][1], pos[p2][2], pos[p1][0], pos[p1][1], pos[p1][2], offset, dist, p1, p2,  potential_lookup[offset + dist], r, step);
        //printf("atom1 %d atom2 %d\n", p1, p2);
        // effectively add zero to energy and force
        if (dist >= number_of_steps) 
            continue; 

        //printf("potential_lookup %lf offset %d dist %d r %lf step %lf number_of_steps %d internal_max %d internal_min %d number_of_types %d \n", potential_lookup[offset + dist], 
          //  offset, dist, r, step, number_of_steps, internal_max, internal_min, number_of_types);

        
        energy += potential_lookup[offset + dist];

        RealOpenMM dEdR = derivative_lookup[offset + dist];
        dEdR = (r > 0) ? (dEdR/r) : 0;
        force[p1] -= delta*dEdR;
        force[p2] += delta*dEdR;
    }

    //auto t6 = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> fp_ms3 = t6 - t5;
    //std::cerr << "execute function " << fp_ms3.count() << "\n";
    //std::cerr << "energy cpu " << energy << std::endl;
    return energy;
}

void ReferenceCalcKBForceKernel::copyParametersToContext(ContextImpl& context, const KBForce& force, const double dist_cutoff) {

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

void ReferenceCalcKBForceKernel::copyTypes(ContextImpl& context, int *atoms, int num_atoms, vector<vector<int>> bonded_exclusions_matrix, const KBForce& force) {
    return;
}

void calculate_non_bonded_forces() {
    return;
}