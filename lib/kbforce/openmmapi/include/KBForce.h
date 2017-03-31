#ifndef OPENMM_KBFORCE_H_
#define OPENMM_KBFORCE_H_

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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportKB.h"

namespace KBPlugin {

/**
 * This class implements an anharmonic bond force of the form E(r)=k*(r-length)^4.  It exists to
 * serve as an example of how to write plugins.
 */

class KBForce : public OpenMM::Force {
public:
    /**
     * Set global parameters
     */
     void setStep(double step) { this->step = step; }
    /**
     * Get global parameters
     */
     double getStep() const { return this->step; }
    /**
     * Get the number of bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }
    /**
     * Add a bond term to the force.
     *
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the force constant for the bond, measured in kJ/mol/nm^4
     * @return the index of the bond that was added
     */
    //~ int addBond(int particle1, int particle2, double length, double k);
    int addBond(int particle1, int particle2, std::vector<double> &potential, std::vector<double> &derivative);
    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to get parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond, measured in kJ/mol/nm^4
     */
    void getBondParameters(int index, int& particle1, int& particle2, std::vector<double>* &potential, std::vector<double>* &derivative) const;
    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to set parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond, measured in kJ/mol/nm^4
     */
    //~ void setBondParameters(int index, int particle1, int particle2, double length, double k);
    void setBondParameters(int index, int particle1, int particle2, std::vector<double> &potential, std::vector<double> &derivative);
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     */
    void updateParametersInContext(OpenMM::Context& context);
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class BondInfo;
    std::vector<BondInfo> bonds;
    double step;
};

class KBForce::BondInfo {
public:
    int particle1, particle2;
	std::vector<double> *potential, *derivative;
    BondInfo() {
        particle1 = particle2 = -1;
		potential = derivative = 0;
    }
    BondInfo(int particle1, int particle2, std::vector<double> &potential, std::vector<double> &derivative) :
        particle1(particle1), particle2(particle2), potential(&potential), derivative(&derivative) {
    }
};

} // namespace KBPlugin

#endif /*OPENMM_KBFORCE_H_*/
