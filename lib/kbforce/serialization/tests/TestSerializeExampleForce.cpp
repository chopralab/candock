/* -------------------------------------------------------------------------- *
 *                                OpenMMExample                                 *
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

#include "ExampleForce.h"
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

extern "C" void registerExampleSerializationProxies();

void testSerialization() {
    // Create a Force.

    ExampleForce force;
    vector<double> kbPotential,
					kbDeriv;
	kbPotential.push_back(50);
	kbPotential.push_back(5);
	kbPotential.push_back(-1);
	kbPotential.push_back(-0.5);
	kbDeriv.push_back(60); 
	kbDeriv.push_back(-20);
	kbDeriv.push_back(0);  
	kbDeriv.push_back(2.0);
    //~ force.addBond(0, 1, 1.0, 2.0);
    //~ force.addBond(0, 2, 2.0, 2.1);
    //~ force.addBond(2, 3, 3.0, 2.2);
    //~ force.addBond(5, 1, 4.0, 2.3);
    force.addBond(0, 1, kbPotential, kbDeriv);
    force.addBond(0, 2, kbPotential, kbDeriv);
    force.addBond(2, 3, kbPotential, kbDeriv);
    force.addBond(5, 1, kbPotential, kbDeriv);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<ExampleForce>(&force, "Force", buffer);
    ExampleForce* copy = XmlSerializer::deserialize<ExampleForce>(buffer);

    // Compare the two forces to see if they are identical.

    ExampleForce& force2 = *copy;
    ASSERT_EQUAL(force.getNumBonds(), force2.getNumBonds());
    for (int i = 0; i < force.getNumBonds(); i++) {
        int a1, a2, b1, b2;
        //~ double da, db, ka, kb;
        vector<double> *pa, *pb, *da, *db;
        //~ force.getBondParameters(i, a1, a2, da, ka);
        //~ force2.getBondParameters(i, b1, b2, db, kb);
        force.getBondParameters(i, a1, a2, pa, da);
        force2.getBondParameters(i, b1, b2, pb, db);
        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        //~ ASSERT_EQUAL(da, db);
        //~ ASSERT_EQUAL(ka, kb);
        ASSERT_EQUAL(pa, pb);
        ASSERT_EQUAL(da, db);
    }
}

int main() {
    try {
        registerExampleSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
