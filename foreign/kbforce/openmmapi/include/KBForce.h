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
#include <set>
#include "internal/windowsExportKB.h"

namespace KBPlugin {

/**
 * This class implements an anharmonic bond force of the form E(r)=k*(r-length)^4.  It exists to
 * serve as an example of how to write plugins.
 */

class OPENMM_EXPORT_KB KBForce : public OpenMM::Force {
public:
    /**
     * Create an KBForce.
     * 
     * @param number_of_types number of types used to create the knowledge based scoring function
     * @param step            step sized used in knowledge based scoring function
     * @param cutoff          cutoff used in knowledge based scoring function
     */
    KBForce(int number_of_types, double step, double cutoff);
    /**
     * Destroy an KBForce.
     */
    ~KBForce();
    /**
     * Get global parameters
     */
     double getStep() const { return this->step; }
     double getCutoff() const { return cutoff; }
    /**
     * Get the number of bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }
    int getNumTypes() const {
        return number_of_types;
    }
    int getNumInteractions() const {
        return number_of_interactions;
    }
    size_t getNumSteps() const {
        return number_of_steps;
    }

    /**
     * Add bonded data
     * 
     * @param type1       IDATM type of particle1
     * @param type2       IDATM type of particle2
     * @param potential   potential of the bond at increasing distance (by step)
     * @param derivative derivative of the bond at increasing distance (by step)
     */

    void addBondType( int type1, int type2,
                      const std::vector<double> &potential,
                      const std::vector<double> &derivative);

    /**
     * Add a bond term to the force.
     *
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param type1     IDATM type of particle1
     * @param type2     IDATM type of particle2
     * @return the index of the bond that was added
     */
    int addBond(int particle1, int particle2, int type1, int type2);
    
    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to get parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param type1     IDATM type of particle1
     * @param type2     IDATM type of particle2
     */
    void getBondParameters(int index, int& particle1, int& particle2, int& type1, int& type2) const;

    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to set parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param type1     IDATM type of particle1
     * @param type2     IDATM type of particle2
     */

    void setBondParameters(int index, int particle1, int particle2, int type1, int type2);

    /**
     * Remove all bonds
     * 
     */
    void clearBonds() {
            bonds.clear();
    }

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
   
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }

    const double* getPotentialLookup() const {
        return &( potential_lookup_table.front());
    }

    const double* getDerivativeLookup() const {
        return &(derivative_lookup_table.front());
    }

protected:
    OpenMM::ForceImpl* createImpl() const;
private:

    struct BondInfo;
    std::vector<BondInfo> bonds;

    int    number_of_types;
    int    number_of_interactions;
    double step;
    double cutoff;
    size_t number_of_steps;
    
    std::vector<double>  potential_lookup_table;
    std::vector<double> derivative_lookup_table;
    
    std::map<int, int> map_internal_to_idatm;
    std::map<int, int> map_idatm_to_internal;

    void insert( const int a, const int b);

};

struct KBForce::BondInfo {
    int particle1, particle2;
    int type1, type2;

    BondInfo() : particle1(-1), particle2(-1), type1(-1), type2(-1) {}
    BondInfo(int particle1, int particle2, int type1, int type2) :
        particle1(particle1), particle2(particle2), type1(type1), type2(type2) {
    }
};

} // namespace KBPlugin

#endif /*OPENMM_KBFORCE_H_*/
