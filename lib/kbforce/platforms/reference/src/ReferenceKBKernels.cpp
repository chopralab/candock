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

//~ double ReferenceCalcKBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //~ vector<RealVec>& pos = extractPositions(context);
    //~ vector<RealVec>& force = extractForces(context);
    //~ int numBonds = particle1.size();
    //~ double energy = 0;
    //~ 
    //~ // Compute the interactions.
    //~ 
    //~ for (int i = 0; i < numBonds; i++) {
        //~ int p1 = particle1[i];
        //~ int p2 = particle2[i];
        //~ RealVec delta = pos[p1]-pos[p2];
        //~ RealOpenMM r2 = delta.dot(delta);
        //~ RealOpenMM r = sqrt(r2);
        //~ RealOpenMM dr = (r-length[i]);
        //~ RealOpenMM dr2 = dr*dr;
        //~ energy += k[i]*dr2*dr2;
        //~ RealOpenMM dEdR = 4*k[i]*dr2*dr;
        //~ dEdR = (r > 0) ? (dEdR/r) : 0;
        //~ force[p1] -= delta*dEdR;
        //~ force[p2] += delta*dEdR;
    //~ }
    //~ return energy;
//~ }
//~ double ReferenceCalcKBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //~ std::cerr << "starting energy calculation (in execute) " << std::endl;
    //~ vector<RealVec>& pos = extractPositions(context);
    //~ vector<RealVec>& force = extractForces(context);
    //~ int numBonds = particle1.size();
    //~ double energy = 0;
    //~ // Compute the interactions.
    //~ 
    //~ for (int i = 0; i < numBonds; i++) {
        //~ int p1 = particle1[i];
        //~ int p2 = particle2[i];
        //~ RealVec delta = pos[p1]-pos[p2];
        //~ RealOpenMM r2 = delta.dot(delta);
        //~ RealOpenMM r = sqrt(r2);
		//~ int dist = (int) round(r / step);
		//~ energy += (*potential[i])[dist];
 //~ 
		//~ std::cerr << "kb particle1 = " << p1 << " particle2 = " << p2 << " distance = " 
			//~ << r << " index for potential = " << dist << " potential = " << (*potential[i])[dist] 
			//~ << " derivative = " << (*derivative[i])[dist] << std::endl;
 //~ 
        //~ RealOpenMM dEdR = (*derivative[i])[dist];
        //~ dEdR = (r > 0) ? (dEdR/r) : 0;
        //~ force[p1] -= delta*dEdR;
        //~ force[p2] += delta*dEdR;
    //~ }
    //~ std::cerr << "kb energy = " << energy << std::endl;
    //~ return energy;
//~ }

//~ double ReferenceCalcKBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //~ // std::cerr << "starting energy calculation (in execute) " << std::endl;
    //~ dbgmsg("starting energy calculation (in execute) ");
    //~ vector<RealVec>& pos = extractPositions(context);
    //~ vector<RealVec>& force = extractForces(context);
    //~ int numBonds = particle1.size();
    //~ double energy = 0;
//~ 
	//~ dbgmsg("numBonds = " << numBonds);
    //~ // Compute the interactions.
    //~ 
    //~ for (int i = 0; i < numBonds; i++) {
        //~ int p1 = particle1[i];
        //~ int p2 = particle2[i];
        //~ RealVec delta = pos[p1] - pos[p2];
        //~ RealOpenMM r2 = delta.dot(delta);
        //~ // RealOpenMM r = 10 * sqrt(r2);
        //~ RealOpenMM r = sqrt(r2);
		//~ // int dist = (int) round(r / step);
		//~ // int dist = (int) round(r / step);
		//~ int dist = (int) floor(r / step);
		//~ if (dist >= (*potential[i]).size()) continue; // effectively add zero to energy and force
        //~ // RealOpenMM dr = (r-length[i]);
        //~ // RealOpenMM dr2 = dr*dr;
        //~ // energy += k[i]*dr2*dr2;
		//~ // energy += scale * (*potential[i])[dist];
		//~ energy += (*potential[i])[dist];
 //~ 
		//~ // std::cerr << "kb particle1 = " << p1 << " particle2 = " << p2 
			//~ // << " distance = " << r << " index for potential = " << dist 
			//~ // << " potential = " << (*potential[i])[dist] 
			//~ // << " derivative = " << (*derivative[i])[dist] 
			//~ // << " size of potential = " << (*potential[i]).size() << std::endl;
//~ 
		//~ // dbgmsg("kb particle1 = " << p1 << " particle2 = " << p2 
			//~ // << " distance = " << r << " index for potential = " << dist 
			//~ // << " potential = " << (*potential[i])[dist] 
			//~ // << " derivative = " << (*derivative[i])[dist] 
			//~ // << " size of potential = " << (*potential[i]).size());
 //~ 
        //~ // RealOpenMM dEdR = 4*k[i]*dr2*dr;
        //~ RealOpenMM dEdR = scale * (*derivative[i])[dist];
        //~ // RealOpenMM dEdR = (*derivative[i])[dist];
        //~ dEdR = (r > 0) ? (dEdR/r) : 0;
        //~ force[p1] -= delta*dEdR;
        //~ force[p2] += delta*dEdR;
    //~ }
    //~ // std::cerr << "kb energy = " << energy << std::endl;
    //~ dbgmsg("knowledge-based energy is calculated as : energy (" << energy 
		//~ << ") * scale (" << setprecision(20) << scale << ") = " 
		//~ << setprecision(20) << energy * scale);
    //~ // return energy;
    //~ return energy * scale;
//~ }

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
		int dist = (int) floor(r / step);
		if (dist >= (*potential[i]).size()) continue; // effectively add zero to energy and force
		energy += (*potential[i])[dist];
		dbgmsg("kb particle1 = " << p1 << " particle2 = " << p2 
			<< " distance = " << r << " index for potential = " << dist 
			<< " potential = " << (*potential[i])[dist] 
			<< " derivative = " << (*derivative[i])[dist] 
			<< " size of potential = " << (*potential[i]).size());

        RealOpenMM dEdR = (*derivative[i])[dist];
        dEdR = (r > 0) ? (dEdR/r) : 0;
        force[p1] -= delta*dEdR;
        force[p2] += delta*dEdR;
    }
    dbgmsg("knowledge-based energy is calculated as : energy = " << setprecision(20) << energy);
    return energy;
}

void ReferenceCalcKBForceKernel::copyParametersToContext(ContextImpl& context, const KBForce& force) {
    if (force.getNumBonds() != particle1.size())
        throw OpenMMException("updateParametersInContext: The number of KB bonds has changed");
    for (int i = 0; i < force.getNumBonds(); i++) {
        int p1, p2;
        //~ force.getBondParameters(i, p1, p2, length[i], k[i]);
        force.getBondParameters(i, p1, p2, potential[i], derivative[i]);
        if (p1 != particle1[i] || p2 != particle2[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}
