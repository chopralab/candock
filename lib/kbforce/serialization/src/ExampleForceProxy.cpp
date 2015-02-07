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

#include "ExampleForceProxy.h"
#include "ExampleForce.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

ExampleForceProxy::ExampleForceProxy() : SerializationProxy("ExampleForce") {
}

void ExampleForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const ExampleForce& force = *reinterpret_cast<const ExampleForce*>(object);
    SerializationNode& bonds = node.createChildNode("Bonds");
    for (int i = 0; i < force.getNumBonds(); i++) {
        int particle1, particle2;
        //~ double distance, k;
		vector<double> *kbPotential, *kbDeriv;
        //~ force.getBondParameters(i, particle1, particle2, distance, k);
        force.getBondParameters(i, particle1, particle2, kbPotential, kbDeriv);
		/*
		 * for now just do nothing
		 */
        //~ bonds.createChildNode("Bond").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("d", distance).setDoubleProperty("k", k);
    }
}

void* ExampleForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    ExampleForce* force = new ExampleForce();
    try {
        const SerializationNode& bonds = node.getChildNode("Bonds");
        for (int i = 0; i < (int) bonds.getChildren().size(); i++) {
            const SerializationNode& bond = bonds.getChildren()[i];
			/*
			 * for now just do nothing
			 */
            //~ force->addBond(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getDoubleProperty("d"), bond.getDoubleProperty("k"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
