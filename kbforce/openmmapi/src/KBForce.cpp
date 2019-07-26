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

#include "KBForce.h"
#include "internal/KBForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>
#include <algorithm>

using namespace KBPlugin;
using namespace OpenMM;
using namespace std;

KBForce::KBForce(int number_of_types, double step, double cutoff) :
     number_of_types(number_of_types),
     number_of_interactions(number_of_types * (number_of_types + 1) / 2),
     step(step), cutoff(cutoff),
     number_of_steps(static_cast<size_t> (cutoff / step)),
     potential_lookup_table( number_of_steps * number_of_interactions ),
     derivative_lookup_table(number_of_steps * number_of_interactions )
     { }

KBForce::~KBForce() {
}

void KBPlugin::KBForce::addBondType(int type1, int type2, const std::vector<double>& potential, const std::vector<double>& derivative) {
    ASSERT( potential.size() >= number_of_steps);
    ASSERT(derivative.size() >= number_of_steps);

    if ( !map_idatm_to_internal.count(type1) ) {
            insert(type1, map_internal_to_idatm.size());
    }

    if ( !map_idatm_to_internal.count(type2) ) {
            insert(type2, map_internal_to_idatm.size());
    }

    const int internal_min = std::min(map_idatm_to_internal.at(type1), map_idatm_to_internal.at(type2) );
    const int internal_max = std::max(map_idatm_to_internal.at(type1), map_idatm_to_internal.at(type2) );

    ASSERT(internal_min < number_of_types );
    ASSERT(internal_max < number_of_types );

    const int offset = (internal_max + internal_min*(number_of_types-1)-internal_min*(internal_min-1)/2) * number_of_steps;

    std::copy( potential.begin(), potential.begin() + number_of_steps, potential_lookup_table.begin() + offset);
    std::copy(derivative.begin(),derivative.begin() + number_of_steps,derivative_lookup_table.begin() + offset);
}


//Add bonds
int KBPlugin::KBForce::addBond (int particle1, int particle2, int type1, int type2) {
    ASSERT(map_idatm_to_internal.count(type1) != 0 );
    ASSERT(map_idatm_to_internal.count(type2) != 0 );

    //map_idatm_to_internal is a lookup table
    bonds.push_back(BondInfo(particle1, particle2,
                             map_idatm_to_internal.at(type1),
                             map_idatm_to_internal.at(type2)));
    return bonds.size()-1;
}

void KBForce::getBondParameters(int index, int& particle1, int& particle2, int& type1, int& type2) const {
    ASSERT_VALID_INDEX(index, bonds);

    particle1 = bonds[index].particle1;
    particle2 = bonds[index].particle2;
    type1 = bonds[index].type1;
    type2 = bonds[index].type2;
}

void KBPlugin::KBForce::setBondParameters (int index, int particle1, int particle2, int type1, int type2) {
    ASSERT_VALID_INDEX(index, bonds);
    ASSERT(map_idatm_to_internal.count(type1) != 0 );
    ASSERT(map_idatm_to_internal.count(type2) != 0 );

    bonds[index].particle1 = particle1;
    bonds[index].particle2 = particle2;
    bonds[index].type1 = map_idatm_to_internal.at(type1);
    bonds[index].type2 = map_idatm_to_internal.at(type2);
}

ForceImpl* KBForce::createImpl() const {
    return new KBForceImpl(*this);
}

void KBForce::updateParametersInContext(Context& context) {
    dynamic_cast<KBForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), 1.5);
}

void KBForce::copyTypes(Context& context, int *atoms, int num_atoms, vector<vector<int>> bonded_exclusions_matrix) {
    dynamic_cast<KBForceImpl&>(getImplInContext(context)).copyTypes(getContextImpl(context), atoms, num_atoms, bonded_exclusions_matrix);
}

void KBPlugin::KBForce::insert(const int a, const int b) {
    map_idatm_to_internal.insert( std::pair<int,int> (a, b));
    map_internal_to_idatm.insert( std::pair<int,int> (b, a));
}
