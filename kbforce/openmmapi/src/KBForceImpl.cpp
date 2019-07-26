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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "internal/KBForceImpl.h"
#include "KBKernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <iostream>

using namespace KBPlugin;
using namespace OpenMM;
using namespace std;

KBForceImpl::KBForceImpl(const KBForce& owner) : owner(owner) {
	//~ std::cout << "calling constructor of KBForceImpl" << endl;
}

KBForceImpl::~KBForceImpl() {
	//~ std::cout << "calling destructor of KBForceImpl" << endl;
}

void KBForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcKBForceKernel::Name(), context);
    kernel.getAs<CalcKBForceKernel>().initialize(context.getSystem(), owner);
}

double KBForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcKBForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> KBForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcKBForceKernel::Name());
    return names;
}

vector<pair<int, int> > KBForceImpl::getBondedParticles() const {
    int numBonds = owner.getNumBonds();
    vector<pair<int, int> > bonds(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int potential, derivative;
        owner.getBondParameters(i, bonds[i].first, bonds[i].second, potential, derivative);
    }
    return bonds;
}

void KBForceImpl::updateParametersInContext(ContextImpl& context, const double dist_cutoff) {
    kernel.getAs<CalcKBForceKernel>().copyParametersToContext(context, owner, dist_cutoff);
}

void KBForceImpl::copyTypes(ContextImpl& context, int *atoms, int num_atoms, vector<vector<int>> bonded_exclusions_matrix) {
    kernel.getAs<CalcKBForceKernel>().copyTypes(context, atoms, num_atoms, bonded_exclusions_matrix, owner);
}

void KBForceImpl::calculate_non_bonded_forces() {
    kernel.getAs<CalcKBForceKernel>().calculate_non_bonded_forces();
}
