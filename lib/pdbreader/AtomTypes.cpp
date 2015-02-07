#include <vector>
#include <map>
#include <set>

#include "Mol.h"
#include "templates/TAexcept.h"
#include "timeit.h"

namespace molecule {

enum BondOrder { AMBIGUOUS, SINGLE, DOUBLE }; // need SINGLE==1, DOUBLE==2
typedef std::set<Ring> Rings;

static int freeOxygens(std::vector<Atom *> &,
					std::map<Atom *, int> &, bool = false);
static bool aromaticGeometry(const Ring &);
static void makeAssignments(std::set<Bond *> &,
	std::map<Bond *, BondOrder> &, std::map<Bond *, int> &,
	std::vector<std::map<Bond *, int> > *, bool allowCharged=false);
static std::map<Bond *, int> *findBestAssignment(
				std::vector<std::map<Bond *, int> > &,
				std::vector<Ring> &);
static bool isN2plus(std::map<Bond *, int> *, const Atom::Bonds &);
static bool isN3plusOkay(const std::vector<Atom *> &);
static void invertUncertains(std::vector<Atom *> &uncertain,
	std::map<Atom *, Bond *> &uncertain2bond,
	std::map<Bond *, BondOrder> *connected);
static void uncertainAssign(std::vector<Atom *> &uncertain,
	std::map<Atom *, Bond *> &uncertain2bond,
	std::set<Bond *> &bonds, std::map<Bond *, BondOrder> &connected,
	std::vector<std::map<Bond *, int> > *assignments,
	std::vector<std::vector<Atom *> > *assignedUncertains,
	bool allowCharged=false);
static void flipAssign(std::vector<Bond *> &flippable, std::set<Atom *> &atoms,
	std::set<Bond *> &bonds, std::map<Bond *, BondOrder> &connected,
	std::vector<std::map<Bond *, int> > *assignments,
	bool allowCharged=false);


template <class Item>
void
generatePermutations(std::vector<Item *> &items,
				std::vector<std::vector<Item *> > *permutations)
{
	std::vector<typename std::vector<Item *> > lastGen;
	std::vector<typename std::vector<Item *>::iterator> lastRems;
	for (typename std::vector<Item *>::iterator ii = items.begin();
	ii != items.end(); ++ii ) {
		Item *i = *ii;
		std::vector<Item *> itemList;
		itemList.push_back(i);
		lastGen.push_back(itemList);
		lastRems.push_back(ii+1);
	}
	permutations->insert(permutations->end(),
					lastGen.begin(), lastGen.end());

	for (unsigned int i = 2; i <= items.size(); ++i) {
		std::vector<std::vector<Item *> > gen;
		std::vector<typename std::vector<Item *>::iterator> rems;
		typename std::vector<std::vector<Item *> >::iterator gi;
		typename std::vector<typename std::vector<Item *>::iterator>::iterator ri;
		for (gi = lastGen.begin(), ri = lastRems.begin();
		gi != lastGen.end(); ++gi, ++ri) {
			for (typename std::vector<Item *>::iterator ii = *ri;
			ii != items.end(); ++ii) {
				std::vector<Item *> perm = *gi;
				perm.push_back(*ii);
				gen.push_back(perm);
				rems.push_back(ii+1);
			}
		}
		permutations->insert(permutations->end(), gen.begin(),
								gen.end());
		lastGen.swap(gen);
		lastRems.swap(rems);
	}
}

static bool
isN2plus(std::map<Bond *, int> *bestAssignment, const Atom::Bonds &bonds)
{
	int sum = 0, target = 4;
	for (Atom::Bonds::const_iterator bi = bonds.begin(); bi != bonds.end();
	++bi) {
		Bond *b = *bi;
		std::map<Bond *, int>::iterator bai = bestAssignment->find(b);
		if (bai == bestAssignment->end())
			target -= 1;
		else
			sum += bai->second;
	}
	return sum == target;
}

static bool
isOarPlus(std::map<Bond *, int> *bestAssignment, const Atom::Bonds &bonds)
{
	int sum = 0, target = 3;
	for (Atom::Bonds::const_iterator bi = bonds.begin(); bi != bonds.end();
	++bi) {
		Bond *b = *bi;
		std::map<Bond *, int>::iterator bai = bestAssignment->find(b);
		if (bai == bestAssignment->end())
			return false;
		else
			sum += bai->second;
	}
	return sum == target;
}

static bool
isN3plusOkay(const std::vector<Atom *> &primary)
{
	for (std::vector<Atom *>::const_iterator bi = primary.begin();
	bi != primary.end(); ++bi) {
		Atom *bondee = *bi;
		Symbol bondeeType = bondee->idatmType();

		if (bondeeType != "C3" && bondeeType != "H" && bondeeType != "D") {
			return false;
		}
	}
	return true;
}

static void
makeAssignments(std::set<Bond *> &bonds,
	std::map<Bond *, BondOrder> &connected,
	std::map<Bond *,int> &curAssign,
	std::vector<std::map<Bond *, int> > *assignments, bool allowCharged)
{
	Bond *assignTarget = *(bonds.begin());
	bonds.erase(bonds.begin());
	bool assign1okay = true, assign2okay = true;
	// see if this assignment completes the bonds of either connected
	// atom and which assignments work
	const Bond::Atoms &bondAtoms = assignTarget->atoms(); 
	for (Bond::Atoms::const_iterator ai = bondAtoms.begin();
	ai != bondAtoms.end(); ++ai) {
		Atom *end = *ai;
		bool complete = true;
		int sum = 0;
		const Atom::Bonds &atomBonds = end->primaryBonds();
		// implied proton treated the same as ambiguous non-ring bond
		bool hasAmbiguous = atomBonds.size() == 2;
		for (Atom::Bonds::const_iterator bi = atomBonds.begin();
		bi != atomBonds.end(); ++bi) {
			Bond *b = *bi;
			if (b == assignTarget)
				continue;
			if (bonds.find(b) != bonds.end()) {
				complete = false;
				break;
			}
			if (curAssign.find(b) != curAssign.end()) {
				sum += curAssign[b];
				continue;
			}
			if (connected[b] == AMBIGUOUS) {
				sum++;
				hasAmbiguous = true;
			} else {
				sum += connected[b];
			}
			if (connected[b] == DOUBLE)
				assign2okay = false;
		}
		if (atomBonds.size() == 2) {
			if (sum == 2)
				assign2okay = false;
		} else {
			if (sum > 3) {
				assign1okay = assign2okay = false;
				break;
			}
			if (sum == 3)
				assign2okay = false;
		}
		if (!complete)
			continue;
		int element = end->element().number();
		if (element > 20)
			continue;
		int valenceElectrons = (element - 2) % 8;
		if (valenceElectrons >= 4) {
			int chargeMod = 0;
			valenceElectrons += sum;
			if (allowCharged) {
				if (element == 7 && atomBonds.size() == 3) {
					chargeMod = 1;
					hasAmbiguous = true;
				} else if (element == 8 && atomBonds.size() == 2) {
					chargeMod = 1;
					hasAmbiguous = true;
				}
			}
			if (hasAmbiguous) {
				if (valenceElectrons < 6
				|| valenceElectrons > 7 + chargeMod)
					assign1okay = false;
				if (valenceElectrons < 5
				|| valenceElectrons > 6 + chargeMod)
					assign2okay = false;
			} else {
				if (valenceElectrons != 7)
					assign1okay = false;
				if (valenceElectrons != 6)
					assign2okay = false;
			}
		} else {
			valenceElectrons -= sum;
			if (hasAmbiguous) {
				if (valenceElectrons < 1
				|| valenceElectrons > 2)
					assign1okay = false;
				if (valenceElectrons < 2
				|| valenceElectrons > 3)
					assign2okay = false;
			} else {
				if (valenceElectrons != 1)
					assign1okay = false;
				if (valenceElectrons != 2)
					assign2okay = false;
			}
		}
	}
	if (assign1okay) {
		curAssign[assignTarget] = 1;
		if (bonds.size() > 0) {
			makeAssignments(bonds, connected, curAssign,
						assignments, allowCharged);
		} else {
			assignments->push_back(curAssign);
		}
		curAssign.erase(assignTarget);
	}
	if (assign2okay) {
		curAssign[assignTarget] = 2;
		if (bonds.size() > 0) {
			makeAssignments(bonds, connected, curAssign,
						assignments, allowCharged);
		} else {
			assignments->push_back(curAssign);
		}
		curAssign.erase(assignTarget);
	}
	bonds.insert(assignTarget);
}

static std::map<Bond *, int> *
findBestAssignment(std::vector<std::map<Bond *, int> > &assignments,
					std::vector<Ring> &systemRings)
{
	if (assignments.size() == 0)
		return NULL;
	if (assignments.size() == 1)
		return &assignments[0];

	// prefer aromatic if possible (and avoid anti-aromatic)
	std::set<int> okayAssignments;
	std::map<int, int> numPlus;
	int bestAro = 0;
	for (unsigned int i = 0; i < assignments.size(); ++i) {
		std::map<Bond *, int> &assignment = assignments[i];
		std::map<Atom *, int> sumOrders, sumBonds;
		for (std::map<Bond *, int>::iterator ai = assignment.begin();
		ai != assignment.end(); ++ai) {
			Bond *b = (*ai).first;
			int order = (*ai).second;
			const Bond::Atoms &bondAtoms = b->atoms();
			for (Bond::Atoms::const_iterator bi =
			bondAtoms.begin(); bi != bondAtoms.end(); ++bi) {
				Atom *a = *bi;
				sumOrders[a] += order;
				sumBonds[a]++;
			}
		}
		// account for N2+
		for (std::map<Atom *, int>::iterator soi = sumOrders.begin();
		soi != sumOrders.end(); ++soi) {
			Atom *a = (*soi).first;
			if (a->element() == Element::N) {
				if (isN2plus(&assignment, a->primaryBonds()))
					numPlus[i] += 1;
			} else if (a->element() == Element::O) {
				if (isOarPlus(&assignment, a->primaryBonds()))
					numPlus[i] += 1;
			}
		}
		int numAro = 0;
		for (std::vector<Ring>::iterator ri = systemRings.begin();
		ri != systemRings.end(); ++ri) {
			Ring ring = *ri;
			int piElectrons = 0;
			const Ring::Atoms &atoms = ring.atoms();
			for (Ring::Atoms::const_iterator ai = atoms.begin();
			ai != atoms.end(); ++ai) {
				Atom *a = *ai;
				int sum = sumOrders[a];
				int element = a->element().number();
				if (element > 20) {
					// guessing not aromatic
					piElectrons = 0;
					break;
				}
				int valenceElectrons = (element - 2) % 8;
				if (valenceElectrons < 3 || valenceElectrons > 6) {
					// aromatic never possible
					piElectrons = 0;
					break;
				}
				if (valenceElectrons == 4) {
					if (sum == 2 && sumBonds[a] == 2) {
						// other bond is double
						piElectrons = 0;
						break;
					}
					piElectrons++;
				} else if (valenceElectrons > 4) {
					piElectrons += 2 - (sum != 2);
				} else if (a->primaryBonds().size() == 2 && sum == 2) {
					piElectrons++;
				} else {
					// aromatic not possible with this assignment
					piElectrons = 0;
					break;
				}
			}
			if (piElectrons % 4 == 2)
				numAro += atoms.size();
		}
		if (numAro > 0 && numAro >= bestAro) {
			if (numAro > bestAro)
				okayAssignments.clear();
			okayAssignments.insert(i);
			bestAro = numAro;
		}
	}
	if (okayAssignments.size() == 1)
		return &assignments[*okayAssignments.begin()];

	// lowest charge preferred
	std::set<int> nextOkay;
	int lowCharge = -1;
	std::vector<std::map<Bond *, int> >::iterator bestAssignment, ai;
	for (ai=assignments.begin(); ai != assignments.end(); ++ai) {
		int index = ai - assignments.begin();
		if (okayAssignments.size() > 0
		&& okayAssignments.find(index) == okayAssignments.end())
			continue;

		int curCharge = numPlus[index];
		if (lowCharge == -1 || curCharge < lowCharge) {
			nextOkay.clear();
			nextOkay.insert(index);
			lowCharge = curCharge;
		} else if (curCharge == lowCharge)
			nextOkay.insert(index);
	}
	okayAssignments.swap(nextOkay);
	if (okayAssignments.size() == 1)
		return &assignments[*okayAssignments.begin()];

	// evaluate by best fit to bond lengths
	float bestVal = 0.0;
	for (ai=assignments.begin(); ai != assignments.end(); ++ai) {
		if (okayAssignments.size() > 0 && okayAssignments.find(
		ai - assignments.begin()) == okayAssignments.end())
			continue;
		float val = 0.0;
		std::map<Bond *, int> &assignment = *ai;
		int orderSum = 0;
		for (std::map<Bond *, int>::iterator i = assignment.begin();
		i != assignment.end(); ++i) {
			val += (*i).first->sqlength() * (*i).second;
			orderSum += (*i).second;
		}
		val /= orderSum;
		if (bestVal == 0.0 || val < bestVal) {
			bestVal = val;
			bestAssignment = ai;
		}
	}
	return &(*bestAssignment);
}

static int
freeOxygens(std::vector<Atom *> &primary, std::map<Atom *, int> &heavys, bool noHyds)
{
	int freeOxygens = 0;
	for (std::vector<Atom *>::const_iterator bi = primary.begin();
	  bi != primary.end(); ++bi) {
		Atom *bondee = *bi;
		if (bondee->element() == Element::O) {
			if (heavys[bondee] != 1)
				continue;
			if (noHyds && bondee->primaryNeighbors().size() != 1)
				continue;
			freeOxygens++;
		}
	}

	return freeOxygens;
}

static bool
aromaticGeometry(const Ring &r)
{
	// algorithm from:
	//	Crystallographic Studies of Inter- and Intramolecular 
	//	   Interactions Reflected in Aromatic Character of pi-Electron
	//	   Systems
	//	J. Chem. Inf. Comput. Sci., 1993, 33, 70-78
	Real sum = 0.0;
	int bonds = 0;
	for (Ring::Bonds::const_iterator bi = r.bonds().begin();
	bi != r.bonds().end(); ++bi) {
		Bond *b = *bi;
		const Element &e1 = b->atoms()[0]->element();
		const Element &e2 = b->atoms()[1]->element();
		Coord c1 = b->atoms()[0]->coord();
		Coord c2 = b->atoms()[1]->coord();
		Real d = distance(c1, c2), delta;

		if (e1 == Element::C && e2 == Element::C) {
			delta = d - 1.38586;
		} else if ((e1 == Element::C || e2 == Element::C) &&
		  (e1 == Element::N || e2 == Element::N)) {
			delta = d - 1.34148;
		} else
			continue;
		bonds++;
		sum += delta * delta;

	}
	if (bonds == 0)
		return false;
	
	Real homa = 1.0 - (792.0 * sum / bonds);

	if (homa >= 0.5) {
		// aromatic
		return true;
	} else if (bonds * homa < -35.0)
		return false;

	return true;
}

void
computeAtomTypes(Molecule *mol)
{
	TimeIt timeit("computeAtomTypes()");

	// angle values used to discriminate between hybridization states
	const Real angle23val1 = 115.0;
	const Real angle23val1_tmax = 116.5;
	const Real angle23val1_tmin = 113.5;
	const Real angle23val2 = 122.0;
	const Real angle12val = 160.0;

	// bond length cutoffs from hybridization discrimination
	// p3... = pass 3 cutoffs; p4... = pass 4 cutoffs
	const Real p3c1c1 = 1.22 * 1.22;
	const Real p3c2c = 1.41 * 1.41;
	const Real p3c2n = 1.37 * 1.37;
	const Real p3n1c1 = 1.20 * 1.20;
	const Real p3n3c = 1.38 * 1.38;
	const Real p3n3n3 = 1.43 * 1.43;
	const Real p3n3n2 = 1.41 * 1.41;
	const Real p3n1o1 = 1.21 * 1.21;
	const Real p3c1o1 = 1.17 * 1.17;
	const Real p3o2c2 = 1.30 * 1.30;
	const Real p3o2as = 1.685 * 1.685;
	const Real p3o2o3 = 1.338 * 1.338;
	const Real p3s2c2 = 1.76 * 1.76;
	const Real p3s2as = 2.11 * 2.11;
	const Real p4c3c = 1.53 * 1.53;
	const Real p4c3n = 1.46 * 1.46;
	const Real p4c3o = 1.44 * 1.44;
	const Real p4n2c = 1.38 * 1.38;
	const Real p4n2n = 1.32 * 1.32;
	const Real p4c2c = 1.42 * 1.42;
	const Real p4c2n = 1.41 * 1.41;
	const Real p4ccnd = 1.45 * 1.45;

	const Real p7cn2nh = 1.3629;
	const Real p7nn2nh = 1.3337;
	const Real p7on2nh = 1.3485;

	// algorithm based on E.C. Meng / R.A. Lewis paper 
	// "Determination of Molecular Topology and Atomic Hybridization
	// States from Heavy Atom Coordinates", J. Comp. Chem., v12#7, 891-898
	// and on example code from idatm.f implementation by E.C. Meng

	// differences: No boron types.  Double-bonded Npls are split off
	//   as N2.  Sox split into Sxd (sulfoxide), and Son (sulfone).
	//   Carbons in aromatic rings are type Car.  Aromatic oxygens are Oar/Oar+.
	//   Negatively charged oxygens are O2- (planar) and O3- (tetrahedral)
	//   instead of O-.  Sp nitrogens bonded to two atoms are N1+.
	
	// since primaryAtoms() returns a copy...
	Atoms pAtoms = mol->primaryAtoms();

	const Atom::IdatmInfoMap &infoMap = Atom::getIdatmInfoMap();

	// initialize idatm type in Atoms
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;
		a->setComputedIdatmType(a->element().name());
	}

	// if molecule is diamond/nanotube, skip atom typing since the
	// ring finding will take forever
	size_t numBonds = mol->numBonds();
	size_t numAtoms = mol->numAtoms();
	if (numBonds - numAtoms > 100 && numBonds / (float) numAtoms > 1.25)
		return;

	std::map<Atom *, int> heavys; // number of heavy atoms bonded
	std::vector<Atom *> primary; // since primaryNeighbors() returns a copy
	size_t hassigned = 0;

	// "pass 1":  type hydrogens / deuteriums and compute number of
	// heavy atoms connected to each atom
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;
		const Element &element = a->element();
		primary = a->primaryNeighbors();

		if (element.number() == 1) {
			// sort out if it's a hydrogen or deuterium
			bool isHyd = true;
			for (const char *c = a->name().c_str(); *c != '\0';
			  ++c) {
				if (isalpha(*c)) {
					if (*c == 'd' || *c == 'D') {
						isHyd = false;
					}
					break;
				}
			}

			bool bondedToCarbon = false;
			for (std::vector<Atom *>::const_iterator bi =
			  primary.begin(); bi != primary.end(); ++bi) {
			  	Atom *bondee = *bi;
				if (bondee->element() == Element::C) {
					bondedToCarbon = true;
					break;
				}
			}
			  
			
			a->setComputedIdatmType(bondedToCarbon ? (isHyd ?
			  "HC" : "DC") : (isHyd ? "H" : "D"));
			hassigned += 1;
		}

		int heavyCount = 0;
		for (std::vector<Atom *>::const_iterator bi = primary.begin();
		  bi != primary.end(); ++bi) {
			Atom *bondee = *bi;
			if (bondee->element().number() > 1) {
				heavyCount++;
			}
		}
		heavys[a] = heavyCount;
	}

	// "pass 1.5": use templates for "infallible" typing of standard
	// residue types
	std::map<const Atom *, bool> mapped;
	const Molecule::Residues &res = mol->residues();
	for (Molecule::Residues::const_iterator ri = res.begin();
	     ri != res.end(); ++ri) {
		Residue *r = *ri;
		try {
			std::vector<Atom *> templated = r->templateAssign(
			 &Atom::setComputedIdatmType, "idatm", "templates", "idatmres");
			for (std::vector<Atom *>::iterator ai =
			templated.begin(); ai != templated.end(); ++ai)
				mapped[*ai] = true;
		} catch (TA_NoTemplate) {
			// don't care
		} catch (...) {
			throw;
		}
	}

	/*
	std::cerr << "  computeIdatmTypes() hydrogens " << hassigned;
	std::cerr << ", template atoms " << mapped.size();
	std::cerr << ", unassigned " << numAtoms - (hassigned + mapped.size()) << std::endl;
	*/

	if (hassigned + mapped.size() == numAtoms)
		return;		// All atoms assigned.

	// "pass 2": elements that are typed only by element type
	// and valences > 1
	std::map<Atom *, int> redo;
	std::set<Atom *> ambiguousVal2Cs;
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;
		if (mapped[a])
			continue;
		const Element &element = a->element();

		// undifferentiated types
		if (element >= Element::He && element <= Element::Be
		  || element >= Element::Ne && element <= Element::Si
		  || element >= Element::Cl) {
		  	a->setComputedIdatmType(element.name());
			continue;
		}

		// valence 4
		//	C must be sp3 (C3)
		//	N must be part of a quaternary amine (N3+) 
		//	P must be part of a phosphate (Pac), a P-oxide (Pox)
		//		or a quaternary phosphine (P3+)
		//	S must be part of a sulfate, sulfonate or sulfamate
		//		(Sac), or sulfone (Son)
		primary = a->primaryNeighbors();
		if (primary.size() == 4) {
			if (element == Element::C) {
				a->setComputedIdatmType("C3");
			} else if (element == Element::N) {
				a->setComputedIdatmType("N3+");
			} else if (element == Element::P) {
				int freeOxys = freeOxygens(primary, heavys);
				if (freeOxys >= 2)
					a->setComputedIdatmType("Pac");
				else if (freeOxys == 1)
					a->setComputedIdatmType("Pox");
				else
					a->setComputedIdatmType("P3+");
			} else if (element == Element::S) {
				int freeOxys = freeOxygens(primary, heavys);
				if (freeOxys >= 3) {
					a->setComputedIdatmType("Sac");
				} else if (freeOxys >= 1) {
					a->setComputedIdatmType("Son");
				} else {
					a->setComputedIdatmType("S");
				}
			}
		}

		// valence 3
		// calculate the three bond angles and average them;
		// since hydrogens may be missing, cannot count on valence
		// to determine the hybridization state.  Average bond angle
		// assists in discriminating hybridization
		//	C may be sp3 (C3), sp2 (C2), or part of a carboxylate
		//		(Cac)
		//	N may be sp3 (N3), sp2, or planar (as in amides and
		//		aniline deriviatives), or part of a nitro
		//		group (Ntr)
		//	S may be, depending on oxidation state, sulfoxide (Sxd)
		//		or S3+
		else if (primary.size() == 3) {
			Real avgAngle = 0.0;
			for (int n1 = 0; n1 < 3; ++n1) {
				for (int n2 = n1 + 1; n2 < 3; ++n2) {
					avgAngle += angle(
					  primary[n1]->coord(), a->coord(),
					  primary[n2]->coord());
				}
			}
			avgAngle /= 3.0;

			if (element == Element::C) {
				bool c3 = false;
				if (avgAngle < angle23val1_tmin)
					c3 = true;
				else if (avgAngle < angle23val1_tmax) {
					Real minSqDist = -1.0;
					for (int n1 = 0; n1 < 3; ++n1) {
						Real sqd = a->coord().sqdistance(primary[n1]->coord());
						if (minSqDist < 0.0 || sqd < minSqDist)
							minSqDist = sqd;
					}
					if (minSqDist > p4c3c)
						c3 = true;
					else if (minSqDist > p4c2c && avgAngle < angle23val1)
						c3 = true;
				}
				if (c3)
					a->setComputedIdatmType("C3");
				else
					a->setComputedIdatmType(freeOxygens(
						primary, heavys, true) >= 2
						? "Cac" : "C2");
			} else if (element == Element::N) {
				if (avgAngle < angle23val1)
					a->setComputedIdatmType("N3");
				else
					a->setComputedIdatmType(freeOxygens(
					  primary, heavys) >= 2 ? "Ntr":"Npl");
			} else if (element == Element::S) {
				bool hasOxy = false;
				for (int i = 0; i < 3; ++i) {
					if (primary[i]->element() ==
					  Element::O) {
					  	hasOxy = true;
						break;
					}
				}
				a->setComputedIdatmType(hasOxy ? "Sxd" : "S3+");
			}
		}

		// valence 2
		// calculate the bond angle and assign a tentative atom
		// type accordingly (a single angle is often not a good
		// indicator of type).  Mark these atoms for further
		// analysis by putting a non-zero value for them in the
		// 'redo' array.
		//	C may be sp3 (C3), sp2 (C2), or sp (C1)
		//	N may be sp3 (N3), sp2 or planar (Npl), or sp (N1+)
		//	O and S are sp3 (O3 and S3, respectively)
		else if (primary.size() == 2) {
			Point coord[2];
			int coordInd = 0;
			for (std::vector<Atom *>::const_iterator bi =
			  primary.begin(); bi != primary.end(); ++bi) {
			  	Atom *other = *bi;
				coord[coordInd++] = other->coord();
			}
			Real ang = angle(coord[0], a->coord(), coord[1]);

			if (element == Element::C) {
				if (ang < angle23val1) {
					a->setComputedIdatmType("C3");
					redo[a] = 1;
					if (ang > angle23val1_tmin)
						ambiguousVal2Cs.insert(a);
				} else if (ang < angle12val) {
					a->setComputedIdatmType("C2");
					if (ang < angle23val2) {
						redo[a] = 3;
					} else {
						// allow ring bond-order code
						// to change this assignment
						redo[a] = -1;
					}
					if (ang < angle23val1_tmax)
						ambiguousVal2Cs.insert(a);
				} else {
					a->setComputedIdatmType("C1");
				}
			} else if (element == Element::N) {
				if (ang < angle23val1) {
					a->setComputedIdatmType("N3");
					redo[a] = 2;
				} else {
					a->setComputedIdatmType(
					  ang < angle12val ?  "Npl" : "N1+");
				}
			} else if (element == Element::O) {
				a->setComputedIdatmType("O3");
			} else if (element == Element::S) {
				a->setComputedIdatmType("S3");
			}
		}
	}

	// "pass 3": determine types of valence 1 atoms.  These were typed
	// by element only in previous pass, but can be typed more accurately
	// now that the atoms they are bonded to have been typed.  Bond
	// lengths are used in this pass.  
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		primary = a->primaryNeighbors();
		if (primary.size() != 1)
			continue;
		
		Atom *bondee = *(primary.begin());
		Real sqlen = sqdistance(bondee->coord(), a->coord());
		Symbol bondeeType = bondee->idatmType();

		
		if (a->idatmType() == "C") {
			if (mapped[a])
				continue;
			if (sqlen <= p3c1c1 && bondeeType == "C1") {
				a->setComputedIdatmType("C1");
			} else if (sqlen <= p3c2c &&
			  bondee->element() == Element::C) {
				a->setComputedIdatmType("C2");
			} else if (sqlen <= p3c2n &&
			  bondee->element() == Element::N) {
				a->setComputedIdatmType("C2");
			} else if (sqlen <= p3c1o1 &&
			  bondee->element() == Element::O &&
			  bondee->primaryNeighbors().size() == 1) {
			  	a->setComputedIdatmType("C1-");
			} else {
				a->setComputedIdatmType("C3");
			}
		} else if (a->idatmType() == "N") {
			if (mapped[a])
				continue;
			if ((sqlen <= p3n1c1 && bondeeType == "C1" ||
			  bondeeType == "N1+") || (sqlen < p3n1o1 &&
			  bondee->element() == Element::O)) {
				a->setComputedIdatmType("N1");
			} else if (sqlen > p3n3c &&
			  (bondeeType == "C2" || bondeeType == "C3")) {
				a->setComputedIdatmType("N3");
			} else if ((sqlen > p3n3n3 && bondeeType == "N3") ||
			  (sqlen > p3n3n2 && bondeeType == "Npl")) {
				a->setComputedIdatmType("N3");
			} else if (bondee->element() == Element::C ||
			  bondee->element() == Element::N) {
				a->setComputedIdatmType("Npl");
			} else {
				a->setComputedIdatmType("N3");
			}
		} else if (a->idatmType() == "O") {
			if (bondeeType == "Cac" || bondeeType == "Ntr" ||
			  bondeeType == "N1+") {
				if (!mapped[a])
					a->setComputedIdatmType("O2-");
			} else if (bondeeType == "Pac" || bondeeType == "Sac"
			  || bondeeType == "N3+" || bondeeType == "Pox"
			  || bondeeType == "Son" || bondeeType == "Sxd") {
				if (mapped[a])
					continue;
				a->setComputedIdatmType("O3-");

				// pKa of 3rd phosphate oxygen is 7...
				if (bondeeType != "Pac")
					continue;
				std::vector<Atom *> bondeePrimary =
						bondee->primaryNeighbors();
				int oxys = 0;
				for (std::vector<Atom *>::const_iterator bpi =
				bondeePrimary.begin(); bpi !=
				bondeePrimary.end(); ++bpi) {
					Atom *bp = *bpi;
					if (bp->element() == Element::O
					&& bp->primaryNeighbors().size() == 1)
						oxys += 1;
				}
				if (oxys < 3)
					continue;
				// if this bond is 0.05 A longer than
				// the other P-O bonds, assume OH
				Real len = distance(bondee->coord(),
								a->coord());
				bool longer = true;
				for (std::vector<Atom *>::const_iterator bpi =
				bondeePrimary.begin(); bpi !=
				bondeePrimary.end(); ++bpi) {
					Atom *bp = *bpi;
					if (bp == a
					|| bp->primaryNeighbors().size() > 1
					|| bp->element() != Element::O)
						continue;
					if (len < distance(bondee->coord(),
					bp->coord()) + 0.05) {
						longer = false;
						break;
					}
				}
				if (longer)
					a->setComputedIdatmType("O3");
			} else if (sqlen <= p3c1o1 &&
			  bondee->element() == Element::C &&
			  bondee->primaryNeighbors().size() == 1) {
			  	if (!mapped[a])
					a->setComputedIdatmType("O1+");
			} else if (sqlen <= p3o2c2 &&
			  bondee->element() == Element::C) {
				if (!mapped[a])
					a->setComputedIdatmType("O2");
				if (!mapped[bondee])
					bondee->setComputedIdatmType("C2");
				redo[bondee] = 0;
			} else if (sqlen <= p3o2as &&
			  bondee->element() == Element::As) {
				if (!mapped[a])
					a->setComputedIdatmType("O2");
			} else if (sqlen <= p3o2o3 &&
			  bondee->element() == Element::O &&
			  bondee->primaryNeighbors().size() == 1) {
				// distinguish oxygen molecule from
				// hydrogen peroxide
				if (!mapped[a])
					a->setComputedIdatmType("O2");
			} else if (sqlen <= p3n1o1 &&
			  bondee->element() == Element::N &&
			  bondee->primaryNeighbors().size() == 1) {
			  	if (!mapped[a])
					a->setComputedIdatmType("O1");
			} else {
				if (!mapped[a])
					a->setComputedIdatmType("O3");
			}
		} else if (a->idatmType() == "S") {
			if (bondee->element() == Element::P) {
				if (!mapped[a])
					a->setComputedIdatmType("S3-");
			} else if (bondeeType == "N1+") {
				if (!mapped[a])
					a->setComputedIdatmType("S2");
			} else if (sqlen <= p3s2c2 &&
			  bondee->element() == Element::C) {
				if (!mapped[a])
					a->setComputedIdatmType("S2");
				if (!mapped[bondee])
					bondee->setComputedIdatmType("C2");
				redo[bondee] = 0;
			} else if (sqlen <= p3s2as &&
			  bondee->element() == Element::As) {
				if (!mapped[a])
					a->setComputedIdatmType("S2");
			} else {
				if (!mapped[a])
					a->setComputedIdatmType("S3");
			}
		}
	}

	// "pass 4": re-examine all atoms with non-zero 'redo' values and
	//   retype them if necessary
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		if (mapped[a])
			redo[a] = 0;

		if (redo[a] == 0)
			continue;
		
		bool c3able = false;
		primary = a->primaryNeighbors();
		for (std::vector<Atom *>::const_iterator bi = primary.begin();
		  bi != primary.end(); ++bi) {
			Atom *bondee = *bi;
			const Element &bondeeElement = bondee->element();
			Real sqlen = sqdistance(bondee->coord(),
								a->coord());

			if (redo[a] == 1) {
				if ((sqlen <= p4c2c && bondeeElement ==
				  Element::C) || (sqlen <= p4c2n &&
				  bondeeElement == Element::N)) {
					a->setComputedIdatmType("C2");
					break;
				}
				if ((sqlen > p4c3c && bondeeElement ==
				  Element::C) || (sqlen > p4c3n &&
				  bondeeElement == Element::N) || (sqlen >
				  p4c3o && bondeeElement == Element::O)) {
					a->setComputedIdatmType("C3");
				}
			} else if (redo[a] == 2) {
				if ((sqlen <= p4n2c && bondeeElement ==
				  Element::C) || (sqlen <= p4n2n &&
				  bondeeElement == Element::N)) {
					// explicit hydrogen(s): N2
					if (heavys[a] < 2)
						a->setComputedIdatmType("N2");
					else
						a->setComputedIdatmType("Npl");
					break;
				}
			} else {
				if ((sqlen <= p4c2c && bondeeElement ==
				  Element::C) || (sqlen <= p4c2n &&
				  bondeeElement == Element::N)) {
					a->setComputedIdatmType("C2");
					c3able = false;
					break;
				}
				if ((sqlen > p4c3c && bondeeElement ==
				  Element::C) || (sqlen > p4c3n &&
				  bondeeElement == Element::N) || (sqlen >
				  p4c3o && bondeeElement == Element::O)) {
					c3able = true;
				}
				if (sqlen > p4ccnd &&
				  bondeeElement == Element::C) {
					c3able = true;
				}
			}
		}
		if (c3able)
			a->setComputedIdatmType("C3");
	}

	// "pass 4.5":  this pass is not in the IDATM paper but is a suggested
	//    improvement mentioned on page 897 of the paper:  find aromatic
	//    ring types.  The method is to:
	//
	//	1) Find all intraresidue rings (actually computed before pass 6)
	//	2) Check that all the atoms of the ring are planar types
	//	3) Check bond lengths around the ring; see if they are
	//		consistent with aromatic bond lengths
	std::set<const Residue *> mappedResidues;
	for (std::map<const Atom *, bool>::const_iterator mi = mapped.begin();
			mi != mapped.end(); ++mi) {
		if ((*mi).second)
			mappedResidues.insert((*mi).first->residue());
	}
	std::set<const Atom *> considerMapped;
	for (std::set<const Residue *>::const_iterator mri = mappedResidues.begin();
			mri != mappedResidues.end(); ++mri) {
		const Residue *r = *mri;
		const Residue::Atoms &rAtoms = r->atoms();
		for (Residue::Atoms::const_iterator rai = rAtoms.begin(); rai != rAtoms.end(); ++rai) {
			considerMapped.insert(*rai);
		}
	}
	int ringLimit = 3;
	Rings rs = mol->rings(false, ringLimit, &considerMapped);
	// std::cerr << "  computeIdatmTypes() size 3 rings " << rs.size() << std::endl;
	if (rs.size() < 20) {
		// not something crazy like an averaged structure...
		ringLimit = 6;
		rs = mol->rings(false, ringLimit, &considerMapped);
		// std::cerr << "  computeIdatmTypes() size 6 rings " << rs.size() << std::endl;
		if (rs.size() < 20) {
			// not something crazy like a nanotube...
			ringLimit = 0;
			rs = mol->rings(false, ringLimit, &considerMapped);
		}
	}
	// screen out rings with definite non-planar types
	Rings planarRings;
	for (Rings::const_iterator ri = rs.begin(); ri != rs.end(); ++ri) {
		const Ring &r = *ri;

		if (r.atoms().size() == 3) {
			for (Ring::Atoms::const_iterator ai = r.atoms().begin();
			ai != r.atoms().end(); ++ai) {
				Atom *a = *ai;
				if (a->element() == Element::C)
					a->setComputedIdatmType("C3");
			}
			continue;
		}
		bool planarTypes = true;
		bool allMapped = true;
		bool allPlanar = true;
		int numOxygens = 0;
		for (Ring::Atoms::const_iterator ai = r.atoms().begin();
		ai != r.atoms().end(); ++ai) {
			const Atom *a = *ai;
			Symbol idatmType = a->idatmType();
			if (a->element() == Element::O)
				numOxygens++;

			primary = a->primaryNeighbors();
			if (primary.size() > 3) {
				allPlanar = planarTypes = false;
				break;
			}

			if (idatmType != "C2" && idatmType != "Npl" &&
			  idatmType != "Sar" && idatmType != "O3" &&
			  idatmType != "S3" && idatmType != "N3" &&
			  idatmType != "Oar" && idatmType != "Oar+" && idatmType != "P" &&
			  idatmType != "Car" && idatmType != "N2" && idatmType != "N2+" &&
			  !(idatmType == "C3" && primary.size()==2)) {
				allPlanar = planarTypes = false;
				break;
			} else if (idatmType == "O3" || idatmType == "S3" ||
			  idatmType == "N3" || idatmType == "C3") {
			  	allPlanar = false;
			}

			if (mapped[a]) {
				if (idatmType == "C3" || idatmType == "O3" ||
				  idatmType == "S3" || idatmType == "N3") {
					allPlanar = planarTypes = false;
					break;
				}
			} else {
				allMapped = false;
			}
		}

		if (!planarTypes)
			continue;
		
		if (allMapped)
			continue;
		if (r.atoms().size() == 5 && numOxygens > 1 && numOxygens < 5)
			continue;
		if (allPlanar || aromaticGeometry(r))
			planarRings.insert(r);
	}
	// find ring systems
	std::set<Ring> seenRings;
	std::vector<std::set<Bond *> > fusedBonds;
	std::vector<std::set<Atom *> > fusedAtoms;
	std::vector<std::vector<Ring> > componentRings;
	std::set<Atom *> ringAssignedNs;
	for (Rings::iterator ri = planarRings.begin(); ri != planarRings.end();
								++ri) {
		const Ring &r = *ri;
		if (seenRings.find(r) != seenRings.end())
			continue;
		std::set<Bond *> systemBonds;
		std::set<Atom *> systemAtoms;
		std::vector<Ring> systemRings;
		std::vector<Ring> queue;
		queue.push_back(r);
		seenRings.insert(r);
		while (queue.size() > 0) {
			Ring qr = queue.back();
			queue.pop_back();
			const Ring::Bonds &bonds = qr.bonds();
			const Ring::Atoms &atoms = qr.atoms();
			systemBonds.insert(bonds.begin(), bonds.end());
			systemAtoms.insert(atoms.begin(), atoms.end());
			systemRings.push_back(qr);
			for (Ring::Bonds::const_iterator bi = bonds.begin();
			bi != bonds.end(); ++bi) {
				const Bond *b = *bi;
				std::vector<Ring> bRings = b->minimumRings();
				for (std::vector<Ring>::iterator bri =
				bRings.begin(); bri != bRings.end(); ++bri) {
					Ring br = *bri;
					if (seenRings.find(br) != seenRings.end())
						continue;
					if (planarRings.find(br) == planarRings.end())
						continue;
					queue.push_back(br);
					seenRings.insert(br);
				}
			}
		}
		fusedBonds.push_back(systemBonds);
		fusedAtoms.push_back(systemAtoms);
		componentRings.push_back(systemRings);
	}
		
	for (unsigned int i=0; i < fusedBonds.size(); ++i) {
		std::set<Bond *> &bonds = fusedBonds[i];
		std::set<Atom *> &atoms = fusedAtoms[i];
		std::vector<Ring> &systemRings = componentRings[i];

		if (atoms.size() > 50) {
			// takes too long to do a massive fused-ring system;
			// assume aromatic
			for (std::set<Atom *>::iterator fai =
			atoms.begin(); fai != atoms.end(); ++fai) {
				Atom *fa = *fai;
				if (fa->element() == Element::C) {
					fa->setComputedIdatmType("Car");
				} else if (fa->element() == Element::O) {
					fa->setComputedIdatmType("Oar");
				}
			}
			continue;
		}

		// skip mapped rings
		if (mapped[*atoms.begin()])
			// since rings shouldn't span residues, one atom mapped => all mapped
			continue;

		// find bonds directly connected to rings
		// and try to judge their order
		std::map<Bond *, BondOrder> connected;
		std::set<std::pair<Atom *, Bond*> > ringNeighbors;
		std::vector<std::pair<Bond *, Atom *> > possiblyAmbiguous;
		for (std::set<Atom *>::iterator ai = atoms.begin();
						ai != atoms.end(); ++ai) {
			Atom *a = *ai;
			primary = a->primaryNeighbors();
			for (std::vector<Atom *>::const_iterator ni =
			primary.begin(); ni != primary.end(); ++ni) {
				Atom *n = *ni;
				if (atoms.find(n) != atoms.end())
					continue;
				Bond *nb = (*a->bondsMap().find(n)).second;
				ringNeighbors.insert(std::pair<Atom *, Bond *>(n, nb));

				if (ambiguousVal2Cs.find(n) != ambiguousVal2Cs.end()) {
					connected[nb] = AMBIGUOUS;
					continue;
				}
				Atom::IdatmInfoMap::const_iterator gi =
					infoMap.find(n->idatmType().str());
				if (gi == infoMap.end()) {
					connected[nb] = SINGLE;
					continue;
				}
				int geom = (*gi).second.geometry;
				if (geom != 3) {
					connected[nb] = SINGLE;
					if (geom == 4
					&& n->primaryNeighbors().size() == 1) {
						possiblyAmbiguous.push_back(
							std::pair<Bond *, Atom*>(nb, n));
					}
					continue;
				}
				if (n->element() == Element::N) {
					// aniline can be planar
					connected[nb] = SINGLE;
					continue;
				}
				// look at neighbors (and grandneighbors)
				bool outsideDouble = false;
				bool ambiguous = (redo[n] == -1);
				std::vector<Atom *> nn = n->primaryNeighbors();
				if (nn.size() == 1) {
					if (a->element() == Element::N && n->element() == Element::O)
						connected[nb] = SINGLE;
					else {
						connected[nb] = DOUBLE;
						possiblyAmbiguous.push_back(
							std::pair<Bond *, Atom*>(nb, n));
					}
					continue;
				}
				for (std::vector<Atom *>::iterator n2i =
				nn.begin(); n2i != nn.end(); ++n2i) {
					Atom *n2 = *n2i;
					if (n2 == a)
						continue;
					gi = infoMap.find(
							n2->idatmType().str());
					if (gi == infoMap.end())
						continue;
					int n2geom = (*gi).second.geometry;
					if (n2geom != 3)
						continue;
					bool allSingle = true;
					std::vector<Atom *> nnn =
							n2->primaryNeighbors();
					for (std::vector<Atom *>::iterator n3i
					= nnn.begin(); n3i != nnn.end(); ++n3i){
						Atom *n3 = *n3i;
						if (n3 == n)
							continue;
						gi = infoMap.find(
							n3->idatmType().str());
						if (gi == infoMap.end())
							continue;
						int n3geom = (*gi).second
								.geometry;
						if (n3geom != 3)
							continue;
						ambiguous = true;
						allSingle = false;
					}
					if (allSingle) {
						outsideDouble = true;
						break;
					}
				}
				if (outsideDouble)
					connected[nb] = SINGLE;
				else if (ambiguous)
					connected[nb] = AMBIGUOUS;
				else
					connected[nb] = DOUBLE;
			}
		}
		std::map<Bond *, int> curAssign;
		std::vector<std::map<Bond *, int> > assignments;
		std::vector<std::vector<Atom *> > assignedUncertains;
		std::map<Atom *, Bond *> uncertain2bond;
		makeAssignments(bonds, connected, curAssign, &assignments);
		if (assignments.size() == 0)
			// try a charged ring
			makeAssignments(bonds, connected, curAssign, &assignments, true);
		else {
			// if there are no aromatic assignments for a ring and the ring
			// has a nitrogen/oxygen, append charged assignments
			bool addCharged = false;
			for (std::vector<Ring>::iterator ri = systemRings.begin();
			ri != systemRings.end(); ++ri) {
				Ring ring = *ri;
				bool hasNO = false;
				for (Ring::Atoms::const_iterator rai = ring.atoms().begin();
				rai != ring.atoms().end(); ++rai) {
					Atom *a = *rai;
					if (a->element() == Element::N || a->element() == Element::O) {
						hasNO = true;
						break;
					}
				}
				if (!hasNO)
					continue;
				bool anyAro = false;
				for (std::vector<std::map<Bond *, int> >::iterator ai =
				assignments.begin(); ai != assignments.end(); ++ai) {
					int bondSum = 0;
					for (Ring::Bonds::const_iterator rbi = ring.bonds().begin();
					rbi != ring.bonds().end(); ++rbi) {
						Bond *b = *rbi;
						bondSum += (*ai)[b];
					}
					int size = ring.bonds().size();
					if (bondSum == size + size/2) {
						anyAro = true;
						break;
					}
				}
				if (!anyAro) {
					addCharged = true;
					break;
				}
			}
			if (addCharged) {
				std::vector<std::map<Bond *, int> > prevAssignments = assignments;
				assignments.clear();
				makeAssignments(bonds, connected, curAssign, &assignments, true);
				assignments.insert(assignments.end(), prevAssignments.begin(),
					prevAssignments.end());
			}
		}
		if (assignments.size() == 0) {
			// see if flipping a possibly-ambiguous bond
			// allows an assignment to be made
			std::vector<std::pair<Real, Bond*> > sortable;
			for (std::vector<std::pair<Bond *, Atom*> >::iterator
			si = possiblyAmbiguous.begin();
			si != possiblyAmbiguous.end(); ++si) {
				Bond *b = (*si).first;
				Atom *a = (*si).second;
				Element e = a->element();
				Real certainty;
				if (e == Element::O) {
					certainty = b->sqlength() - p3o2c2;
				} else if (e == Element::S) {
					certainty = b->sqlength() - p3s2c2;
				} else {
					certainty = 0.0;
				}
				if (certainty < 0.0)
					certainty = 0.0 - certainty;
				sortable.push_back(std::pair<Real, Bond *>(
					certainty, b));
			}
			std::sort(sortable.begin(), sortable.end());
			std::vector<Bond *> flippable;
			for (std::vector<std::pair<Real, Bond*> >::iterator si
			= sortable.begin(); si != sortable.end(); ++si) {
				flippable.push_back((*si).second);
			}
			if (flippable.size() > 0) {
				flipAssign(flippable, atoms, bonds, connected, &assignments);
				if (assignments.size() == 0)
					flipAssign(flippable, atoms, bonds,
						connected, &assignments, true);
			}
		}

		if (assignments.size() == 0) {
			// if adjacent carbons were uncertain (i.e. had
			// "redo" values) try changing their type
			std::vector<Atom *> uncertain;
			for (std::set<std::pair<Atom *, Bond *> >::iterator rni =
			ringNeighbors.begin(); rni != ringNeighbors.end(); ++rni) {
				Atom *rna = (*rni).first;
				Bond *rnb = (*rni).second;
				if (redo.find(rna) != redo.end() && redo[rna] != 0) {
					uncertain.push_back(rna);
					uncertain2bond[rna] = rnb;
				}
			}
			if (uncertain.size() > 0) {
				uncertainAssign(uncertain, uncertain2bond, bonds, connected,
							&assignments, &assignedUncertains);
				if (assignments.size() == 0)
					uncertainAssign(uncertain, uncertain2bond, bonds, connected,
								&assignments, &assignedUncertains, true);
			}
		}
		if (assignments.size() == 0) {
			std::cerr << "Cannot find consistent set of"
			" bond orders for ring system containing atom "
			<< (*atoms.begin())->oslIdent() << "\n";
			continue;
		}

		std::map<Bond *, int> *bestAssignment = findBestAssignment(
								assignments, systemRings);
		if (bestAssignment != NULL && assignedUncertains.size() > 0) {
			unsigned int baIndex = std::find(assignments.begin(),
				assignments.end(), *bestAssignment) - assignments.begin();
			invertUncertains(assignedUncertains[baIndex],
							uncertain2bond, &connected);
		}

		// see if individual rings are aromatic -- if not
		// then assign types as per best assignment
		for (std::vector<Ring>::iterator ri = systemRings.begin();
		ri != systemRings.end(); ++ri) {
			Ring ring = *ri;
			int minFreeElectrons = 0, maxFreeElectrons = 0;
			std::vector<Bond *> otherFused;
			const Ring::Bonds &bonds = ring.bonds();
			std::set<Bond *> ringBonds;
			ringBonds.insert(bonds.begin(), bonds.end());
			const Ring::Atoms &atoms = ring.atoms();
			bool aro = true;
			for (Ring::Atoms::const_iterator ai = atoms.begin();
			ai != atoms.end(); ++ai) {
				Atom *a = *ai;
				const Atom::Bonds &a_bonds = a->primaryBonds();
				minFreeElectrons++; // a few exceptions below
				if (a_bonds.size() == 2) {
					int element = a->element().number();
					if (element > 20) {
						maxFreeElectrons += 2;
						continue;
					}
					int valence = (element - 2) % 8;
					if (valence < 3 || valence == 7) {
						aro = false;
						break;
					}
					if (valence < 5)
						maxFreeElectrons++;
					else if (valence == 5) {
						if (bestAssignment == NULL)
							maxFreeElectrons += 2;
						else {
							int sum = (*bestAssignment)[a_bonds[0]] + (*bestAssignment)[a_bonds[1]];
							if (sum == 2) {
								minFreeElectrons++;
								maxFreeElectrons += 2;
							} else
								maxFreeElectrons++;
						}
					} else {
						if (bestAssignment != NULL
						&& (element == 8 && isOarPlus(bestAssignment, a_bonds))) {
								maxFreeElectrons++;
						} else {
							minFreeElectrons++;
							maxFreeElectrons += 2;
						}
					}
				} else if (a_bonds.size() == 3) {
					BondOrder bo = AMBIGUOUS;
					Bond *outBond;
					for (Atom::Bonds::const_iterator abi =
					a_bonds.begin(); abi != a_bonds.end();
					++abi) {
						Bond *b = *abi;
						if (ringBonds.find(b)
						!= ringBonds.end())
							continue;
						outBond = b;
						if (connected.find(b)
						!= connected.end()) {
							bo = connected[b];
							break;
						}
						for (std::vector<std::map<Bond *, int> >::iterator ai = assignments.begin(); ai != assignments.end(); ++ai) {
							int assignment=(*ai)[b];
							if (bo == AMBIGUOUS) {
								bo = (BondOrder) assignment;
							} else if (bo
							!= assignment) {
								bo = AMBIGUOUS;
								break;
							}
						}
						break;
					}
					int element = a->element().number();
					if (element > 20) {
						maxFreeElectrons += 2;
						continue;
					}
					int valence = (element - 2) % 8;
					if (valence < 4 || valence > 5) {
						aro = false;
						break;
					}
					if (valence == 4) {
						if (bo == DOUBLE) {
							if (outBond->otherAtom(a)->primaryNeighbors().size() > 1) {
								bool isFused = false;
								for (std::vector<Ring>::iterator sri = systemRings.begin(); sri != systemRings.end(); ++sri) {
									if (sri == ri) {
										continue;
									}
									Ring sring = *sri;
									const Ring::Atoms &srAtoms = sring.atoms();
									if (srAtoms.find(a) != srAtoms.end()) {
										isFused = true;
										break;
									}
								}
								if (!isFused) {
									aro = false;
									break;
								}
							} else {
								aro = false;
								break;
							}
						}
						maxFreeElectrons += 1;
					} else if (bo == DOUBLE ||
						// ... or N2+ or Oar+
					(bestAssignment != NULL
					&& (element == 7 && isN2plus(bestAssignment, a_bonds))))
						maxFreeElectrons += 1;
					else {
						minFreeElectrons++;
						maxFreeElectrons += 2;
					}
				} else {
					aro = false;
					break;
				}
			}
			int aroEval;
			if (aro) {
				aroEval = (minFreeElectrons-2) % 4;
				if ((aroEval != 0) && (aroEval +
				(maxFreeElectrons - minFreeElectrons) < 4))
					aro = false;
			}
			// assign aromatic types if aro (N2/Npl depends
			// on number of needed electrons)
			std::set<Atom *> protonatableNs;
			for (Ring::Atoms::const_iterator ai = atoms.begin();
			ai != atoms.end(); ++ai) {
				Atom *a = *ai;
				Element e = a->element();
				if (e == Element::N) {
					if (a->primaryBonds().size() == 2)
						protonatableNs.insert(a);
					else if (bestAssignment != NULL) {
						// distinguish N2+/Npl
						if (aro && isN2plus(bestAssignment, a->primaryBonds())) {
							// non-ring bond == 0
							a->setComputedIdatmType("N2+");
							// type bonded isolated O as O3-
							if (a->primaryBonds().size() == 3) {
								Atom::BondsMap::const_iterator end = a->bondsMap().end();
								for (Atom::BondsMap::const_iterator bmi = a->bondsMap().begin();
								bmi != end; ++bmi) {
									Atom *nb = (*bmi).first;
									if (nb->element() == Element::O
									&& nb->bondsMap().size() == 1)
										nb->setComputedIdatmType("O3-");
								}
							}
						}
					}
					continue;
				}
				if (!aro)
					continue;
				if (e == Element::C) {
					a->setComputedIdatmType("Car");
				} else if (e == Element::O) {
					if (isOarPlus(bestAssignment, a->primaryBonds()))
						a->setComputedIdatmType("Oar+");
					else
						a->setComputedIdatmType("Oar");
				} else if (e == Element::S) {
					a->setComputedIdatmType("Sar");
				} else if (e == Element::P) {
					if (a->primaryBonds().size() == 2)
						if (aroEval % 4 != 0)
							aroEval++;
				} else if (e.number() > 20) {
					if (aroEval % 4 != 0)
						aroEval++;
				}
			}
			if (bestAssignment == NULL && aro) {
				// N2/Npl depends on number of needed electrons
				while (aroEval % 4 != 0 && protonatableNs.size() > 0) {
					aroEval++;
					Atom *longestN = NULL;
					float bestVal = 0.0;
					for (std::set<Atom *>::iterator ai =
					protonatableNs.begin(); ai !=
					protonatableNs.end(); ++ai) {
						Atom *a = *ai;
						float val = 0.0;
						const Atom::Bonds &a_bonds = a->primaryBonds();
						for (Atom::Bonds::const_iterator bi =
						a_bonds.begin(); bi != a_bonds.end();
						++bi) {
							Bond *b = *bi;
							val += b->length();
						}
						if (longestN == NULL || val > bestVal) {
							longestN = a;
							bestVal = val;
						}
						// avoid retyping in pass 7
						redo[a] = -7;
					}
					longestN->setComputedIdatmType("Npl");
					protonatableNs.erase(longestN);
					ringAssignedNs.insert(longestN);

				}
				for (std::set<Atom *>::iterator ai =
				protonatableNs.begin(); ai != protonatableNs.end();
				++ai) {
					Atom *a = *ai;
					a->setComputedIdatmType("N2");
					ringAssignedNs.insert(a);
				}
			} else if (bestAssignment != NULL) {
				
				// decide if two-bond nitrogens are N2 or Npl
				for (std::set<Atom *>::iterator ai =
				protonatableNs.begin(); ai !=
				protonatableNs.end(); ++ai) {
					Atom *a = *ai;
					int bondSum = 0;
					const Atom::Bonds &a_bonds = a->primaryBonds();
					for (Atom::Bonds::const_iterator bi =
					a_bonds.begin(); bi != a_bonds.end();
					++bi) {
						Bond *b = *bi;
						bondSum += (*bestAssignment)[b];
					}
					if (bondSum == 2)
						a->setComputedIdatmType("Npl");
					else
						a->setComputedIdatmType("N2");
					ringAssignedNs.insert(a);
				}
			}
		}
	}

	// "pass 5": change isolated sp2 carbons to sp3 since it is 
	//   impossible for an atom to be sp2 hybrizided if all its 
	//   neighbors are sp3 hybridized.  In addition, a carbon atom cannot
	//   be double bonded to a carboxylate carbon, phosphate phosphorus,
	//   sulfate sulfur, sulfone sulfur, sulfoxide sulfur, or sp1 carbon.
	//   Addition not in original idatm: Npl/N2+ also
	//
	//   This has now been streamlined to:  must be bonded to an sp2
	//   atom other than carboxylate carbon. Also, if the sp2 carbon
	//   is valence 3 and a neighbor is valence 1, then "trust" the sp2
	//   assignment and instead change the neighbor to sp2.
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		if (a->idatmType() != "C2")
			continue;

		if (mapped[a])
			continue;
		
		bool c2possible = false;
		primary = a->primaryNeighbors();
		std::vector<Atom *> nbValence1;
		for (std::vector<Atom *>::const_iterator bi = primary.begin();
		  bi != primary.end(); ++bi) {
			Atom *bondee = *bi;
			Symbol bondeeType = bondee->idatmType();

			Atom::IdatmInfoMap::const_iterator i = infoMap.find(
							bondeeType.str());
			if (i == infoMap.end())
				continue;
			if ((*i).second.geometry == 3 && bondeeType != "Cac"
					&& bondeeType != "N2+"
					// Npl with two bonds or less may be N2
					&& !(bondeeType == "Npl"
					&& bondee->primaryNeighbors().size() > 2
					// because Ng+ isn't assigned until next pass
					&& heavys[bondee] > 1)) {
				c2possible = true;
			  	break;
			} else if (bondee->primaryNeighbors().size() == 1)
				nbValence1.push_back(bondee);
		}

		if (!c2possible) {
			if (primary.size() == 3 && nbValence1.size() > 0) {
				Atom *best = NULL;
				float bestRatio;
				const char *bestSp2Type;
				for (std::vector<Atom *>::iterator nbi = nbValence1.begin();
				nbi != nbValence1.end(); ++nbi) {
					Atom *nb = *nbi;
					const char *sp2Type;
					Real test;
					if (nb->element() == Element::C) {
						sp2Type = "C2";
						test = p3c2c;
					} else if (nb->element() == Element::O) {
						sp2Type = "O2";
						test = p3o2c2;
					} else if (nb->element() == Element::N) {
						sp2Type = "N2";
						test = p4n2c;
					} else if (nb->element() == Element::S) {
						sp2Type = "S2";
						test = p3s2c2;
					} else
						continue;
					Real sqlen = sqdistance(nb->coord(), a->coord());
					Real ratio = sqlen / test;
					if (best == NULL || ratio < bestRatio) {
						best = nb;
						bestRatio = ratio;
						bestSp2Type = sp2Type;
					}
				}
				if (best != NULL)
					best->setComputedIdatmType(bestSp2Type);
				else
					a->setComputedIdatmType("C3");
			} else
				a->setComputedIdatmType("C3");
		}
	}

	// "pass 6": 
	//   1) make decisions about the charge states of nitrogens.  If a
	//      nitrogen is bonded to sp3 carbons and/or hydrogens and/or
	//      deuteriums only, assume that it is positively charged (the pKa
	//      of its conjugate acid is probably high enough that the
	//      protonated form predominates at physiological pH).  If an sp2
	//      carbon is bonded to three planar nitrogens, it may be part of
	//      a guanidinium group.  Make the nitrogens positively charged
	//      (Ng+) if guanidine or similar structures can be ruled out (if
	//      'noplus' is false).
	//   2) make carboxyl oxygens negatively charged even if the proton is
	//      present (the pKa of the carboxyl group is probably low enough
	//      that the unprotonated form predominates at physiological pH).
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;
		
		primary = a->primaryNeighbors();
		const Element &element = a->element();
		if (element == Element::N && a->idatmType() != "N3+") {
			if (mapped[a])
				continue;
			if (isN3plusOkay(primary))
				a->setComputedIdatmType("N3+");
			
		} else if (a->idatmType() == "C2") {
			int numNpls = 0;
			for (std::vector<Atom *>::const_iterator bi =
			  primary.begin(); bi != primary.end(); ++bi) {
				Atom *bondee = *bi;

				if ((bondee->idatmType() == "Npl"
				  && !mapped[bondee])
				  || bondee->idatmType() == "Ng+")
					// Ng+ possible through template
					// typing
					numNpls++;
			}

			bool noplus = false;
			if (numNpls >= 2) {
				for (std::vector<Atom *>::const_iterator bi =
				  primary.begin(); bi != primary.end(); ++bi) {
					Atom *bondee = *bi;

					if (bondee->idatmType() != "Npl")
						continue;
					if (mapped[bondee])
						continue;
					
					if (bondee->rings(false,
							ringLimit, &considerMapped).size() > 0) {
						noplus = true;
						break;
					}
					bondee->setComputedIdatmType("Ng+");
				}
			}
			if (noplus) {
				for (std::vector<Atom *>::const_iterator bi =
				  primary.begin(); bi != primary.end(); ++bi) {
					Atom *bondee = *bi;
					if (mapped[bondee])
						continue;
					if (bondee->idatmType() == "Ng+")
						bondee->setComputedIdatmType("Npl");
				}
			}
		} else if (a->idatmType() == "Cac") {
			for (std::vector<Atom *>::const_iterator bi =
			  primary.begin(); bi != primary.end(); ++bi) {
				Atom *bondee = *bi;
				if (mapped[bondee])
					continue;
				if (bondee->element() == Element::O &&
				  heavys[bondee] == 1) {
					bondee->setComputedIdatmType("O2-");
				}
			}
		}
	}

	// "pass 7":  a non-IDATM pass:  split off heavy-atom-valence-2
	//	Npls that have no hydrogens as type N2.
	//	Discrimination criteria is the average bond length of the two 
	//	heavy-atom bonds (shorter implies more double-bond character,
	//	thereby no hydrogen).
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		if (mapped[a])
			continue;

		if (a->idatmType() != "Npl")
			continue;
		
		if (heavys[a] != 2)
			continue;

		if (ringAssignedNs.find(a) != ringAssignedNs.end())
			continue;
		
		primary = a->primaryNeighbors();
		if (primary.size() > 2)
			continue;

		Real threshold = 1.0;
		Real harmLen = 1.0;
		Atom *recipient = NULL;
		Real bratio = 1.0;
		for (std::vector<Atom *>::const_iterator bi = primary.begin();
		bi != primary.end(); ++bi) {
			Atom *bondee = *bi;
			
			Real criteria;
			if (bondee->element() == Element::C) {
				criteria = p7cn2nh;
			} else if (bondee->element() == Element::N) {
				criteria = p7nn2nh;
			} else if (bondee->element() == Element::O) {
				if (bondee->primaryNeighbors().size() > 1)
					continue;
				criteria = p7on2nh;
			} else {
				continue;
			}
			threshold *= criteria;
			Real len = distance(bondee->coord(), a->coord());
			harmLen *= len;
			if (len > criteria)
				continue;
			Real ratio = len / criteria;
			if (ratio > bratio)
				continue;
			if (bondee->element() == Element::N && bondee->primaryBonds().size() > 2)
				continue;
			if (bondee->idatmType() == "Car")
				continue;

			std::vector<Atom *> bondeePrimary = bondee->primaryNeighbors();
			bool doubleOkay = true;
			for (std::vector<Atom *>::const_iterator bi = bondeePrimary.begin();
			  bi != bondeePrimary.end(); ++bi) {
				Atom *grandBondee = *bi;
				if (grandBondee == a)
					continue;
				Symbol gbType = grandBondee->idatmType();

				Atom::IdatmInfoMap::const_iterator i = infoMap.find(gbType.str());
				if (i == infoMap.end())
					continue;
				int geom = (*i).second.geometry;
				if (geom > 1 && geom < 4 && heavys[grandBondee] == 1) {
						doubleOkay = false;
						break;
				}
			}
			if (!doubleOkay)
				continue;
			recipient = bondee;
			bratio = ratio;
		}

		if (harmLen < threshold && recipient != NULL) {
			a->setComputedIdatmType("N2");
			if (recipient->element() == Element::C) {
				recipient->setComputedIdatmType("C2");
			} else if (recipient->element() == Element::N) {
				recipient->setComputedIdatmType("N2");
			} else if (recipient->element() == Element::O) {
				recipient->setComputedIdatmType("O2");
			}
		}
	}

	// "pass 8":  another non-IDATM: change planar nitrogens bonded only
	//  SP3 atoms to N3 or N3+.  Change Npls/N3s bonded to sp2 atoms that
	//  are not in turn bonded to sp2 atoms (implying Npl doubled bonded)
	//  to N2, otherwise to Npl.  
	std::map<Atom *, std::vector<Atom *> > bondedSp2s;
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		if (mapped[a])
			continue;

		if (ringAssignedNs.find(a) != ringAssignedNs.end())
			continue;
		
		Symbol idatmType = a->idatmType();
		if (idatmType != "Npl" && idatmType != "N2" && idatmType!="N3")
			continue;
		
		bool aroRing = false;
		const Atom::Rings &aRings = a->rings(false, ringLimit, &considerMapped);
		for (Atom::Rings::const_iterator ri = aRings.begin();
		ri != aRings.end(); ++ ri) {
			if ((*ri)->aromatic()) {
				aroRing = true;
				break;
			}
		}
		if (aroRing)
			continue;

		// any sp2?
		primary = a->primaryNeighbors();
		if (idatmType == "Npl" && primary.size() != 2)
			continue;

		std::vector<Atom *> bsp2list;
		for (std::vector<Atom *>::const_iterator bi = primary.begin();
		bi != primary.end(); ++bi) {
			Atom *bondee = *bi;
			Symbol idatmType = bondee->idatmType();

			aroRing = false;
			const Atom::Rings &bRings = bondee->rings(false,
								ringLimit, &considerMapped);
			for (Atom::Rings::const_iterator ri = bRings.begin();
			ri != bRings.end(); ++ ri) {
				if ((*ri)->aromatic()) {
					aroRing = true;
					break;
				}
			}
			if (aroRing) {
				if (heavys[a] == 1) { // aniline
					a->setComputedIdatmType("Npl");
					break;
				}
				continue;
			}

			Atom::IdatmInfoMap::const_iterator i = infoMap.find(
							idatmType.str());
			if (i == infoMap.end() || (*i).second.geometry != 3
					|| bondee->idatmType() == "Npl")
				continue;
			bsp2list.push_back(bondee);
		}
		bondedSp2s[a] = bsp2list;
	}

	// order typing by easiest-figure-out (1 sp2 bonded) to hardest (3)
	// good test cases:  1CY in 3UM8; WRA in 1J3I
	for (unsigned int i=1; i<4; ++i) {
		for (std::map<Atom *, std::vector<Atom *> >::const_iterator sp2i = bondedSp2s.begin();
		sp2i != bondedSp2s.end(); ++sp2i) {
			const std::vector<Atom *> &sp2s = (*sp2i).second;
			if (sp2s.size() != i)
				continue;
			Atom *a = (*sp2i).first;
			bool anySP2 = false;

			for (std::vector<Atom *>::const_iterator bi = sp2s.begin();
			bi != sp2s.end(); ++bi) {
				Atom *bondee = *bi;

				anySP2 = true;
				bool remoteSP2 = false;
				std::vector<Atom *> grandPrimary = bondee->primaryNeighbors();
				for (std::vector<Atom *>::const_iterator gbi =
				grandPrimary.begin(); gbi != grandPrimary.end(); ++gbi){
					Atom *grandBondee = *gbi;
					if (grandBondee == a)
						continue;
					Atom::IdatmInfoMap::const_iterator gi =
					infoMap.find(grandBondee->idatmType().str());
					if (gi == infoMap.end()
					|| (*gi).second.geometry == 3) {
						if (grandBondee->idatmType() != "Car"
						&& grandBondee->idatmType() != "Npl"
						) {
						//&& !(grandBondee->idatmType() == "Npl"
						//& grandBondee->primaryBonds().size()
						//= 2)) {
							remoteSP2 = true;
							break;
						}
					}
				}
				if (!remoteSP2) {
					a->setComputedIdatmType("N2");
					break;
				}
				// a remote sp2 atom doesn't necessarily mean Npl
				// (see N1 in MX1 of 2aio), so no else clause
			}
			if (!anySP2) {
				int hvys = heavys[a];
				if (hvys > 1)
					a->setComputedIdatmType("N3");
				else if (hvys == 1)
					a->setComputedIdatmType(isN3plusOkay(primary) ? "N3+" : "N3");
				else
					a->setComputedIdatmType("N3+");
			}
		}
	}

	// "pass 9":  another non-IDATM pass and analogous to pass 8:
	//  change O3 bonded only to non-Npl sp2 atom not in turn bonded
	//  to non-Npl sp2 to O2.
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		if (mapped[a])
			continue;

		Symbol idatmType = a->idatmType();
		if (idatmType != "O3")
			continue;
		
		primary = a->primaryNeighbors();
		if (primary.size() != 1)
			continue;

		// any sp2?
		Atom *bondee = *(primary.begin());
		Symbol bondeeType = bondee->idatmType();

		bool aroRing = false;
		const Atom::Rings &bRings = bondee->rings(false, ringLimit, &considerMapped);
		for (Atom::Rings::const_iterator ri = bRings.begin();
		ri != bRings.end(); ++ ri) {
			if ((*ri)->aromatic()) {
				aroRing = true;
				break;
			}
		}
		if (aroRing) {
			// can't be O2
			continue;
		}

		Atom::IdatmInfoMap::const_iterator i = infoMap.find(
							bondeeType.str());
		if (i == infoMap.end() || (*i).second.geometry != 3)
			continue;
		bool remoteSP2 = false;
		std::vector<Atom *> grandPrimary = bondee->primaryNeighbors();
		for (std::vector<Atom *>::const_iterator gbi =
		grandPrimary.begin(); gbi != grandPrimary.end(); ++gbi) {
			Atom *grandBondee = *gbi;
			if (grandBondee == a)
				continue;
			Atom::IdatmInfoMap::const_iterator gi =
				infoMap.find(grandBondee->idatmType().str());
			if (gi == infoMap.end() || (*gi).second.geometry == 3) {
				if (grandBondee->idatmType() != "Car"
				&& grandBondee->idatmType() != "Npl") {
					remoteSP2 = true;
					break;
				}
			}
		}
		if (!remoteSP2)
			a->setComputedIdatmType("O2");
	}

	// "pass 10":  another non-IDATM pass. Ensure nitrate ions are N2/O2-
	for (Atoms::iterator ai = pAtoms.begin(); ai != pAtoms.end(); ++ai) {
		Atom *a = *ai;

		if (mapped[a])
			continue;

		if (a->element() != Element::N)
			continue;
		
		primary = a->primaryNeighbors();
		if (primary.size() != 2)
			continue;

		bool bondersOkay = true;
		for (Atoms::const_iterator pi = primary.begin(); pi != primary.end(); ++pi) {
			Atom *bondee = *pi;
			
			if (bondee->element() != Element::O
			|| bondee->primaryNeighbors().size() != 1) {
				bondersOkay = false;
				break;
			}
		}

		if (bondersOkay) {
			a->setComputedIdatmType("N2");
			for (Atoms::const_iterator pi = primary.begin(); pi != primary.end(); ++pi)
				(*pi)->setComputedIdatmType("O2-");
		}
	}
}

static void
invertUncertains(std::vector<Atom *> &uncertain,
	std::map<Atom *, Bond *> &uncertain2bond,
	std::map<Bond *, BondOrder> *connected)
{
	for (std::vector<Atom *>::iterator ui = uncertain.begin();
	ui != uncertain.end(); ++ui) {
		Atom *a = *ui;
		Bond *b = uncertain2bond[a];
		(*connected)[b] = (BondOrder) (3 - (*connected)[b]);
		if (a->idatmType() == "C3")
			a->setComputedIdatmType("C2");
		else if (a->idatmType() == "C2")
			a->setComputedIdatmType("C3");
		else if (a->idatmType() == "Npl")
			a->setComputedIdatmType("N2");
		else if (a->idatmType() == "N2")
			a->setComputedIdatmType("Npl");
		else
			std::cerr << "Unknown redo atom type: "
						<< a->idatmType() << "\n";
	}
}

static void
uncertainAssign(std::vector<Atom *> &uncertain,
	std::map<Atom *, Bond *> &uncertain2bond,
	std::set<Bond *> &bonds, std::map<Bond *, BondOrder> &connected,
	std::vector<std::map<Bond *, int> > *assignments,
	std::vector<std::vector<Atom *> > *assignedUncertains,
	bool allowCharged)
{
	std::map<Bond *, int> curAssign;
	std::vector<std::vector<Atom *> > permutations;
	generatePermutations<Atom>(uncertain, &permutations);
	// if we find an assignment involving changing N atoms, have to
	// try all permutations that involve changing no more than N atoms
	unsigned int limit = uncertain.size();
	for (std::vector<std::vector<Atom *> >::iterator pi =
	permutations.begin(); pi != permutations.end(); ++pi) {
		std::vector<Atom *> &uncertain = *pi;
		if (uncertain.size() > limit)
			break;
		invertUncertains(uncertain, uncertain2bond, &connected);
		unsigned int numPrevAssigned = assignments->size();
		makeAssignments(bonds, connected, curAssign, assignments, allowCharged);
		int increase = assignments->size() - numPrevAssigned;
		if (increase) {
			limit = uncertain.size();
			for (int i=0; i < increase; ++i) {
				assignedUncertains->push_back(uncertain);
			}
		}
		invertUncertains(uncertain, uncertain2bond, &connected);
	}
}

static void
flipAssign(std::vector<Bond *> &flippable, std::set<Atom *> &atoms,
	std::set<Bond *> &bonds, std::map<Bond *, BondOrder> &connected,
	std::vector<std::map<Bond *, int> > *assignments, bool allowCharged)
{
	std::map<Bond *, int> curAssign;
	std::vector<std::vector<Bond *> > permutations;
	generatePermutations<Bond>(flippable, &permutations);
	for (std::vector<std::vector<Bond *> >::iterator pi =
	permutations.begin(); pi != permutations.end(); ++pi) {
		std::vector<Bond *> &flip = *pi;
		for (std::vector<Bond *>::iterator bi =
		flip.begin(); bi != flip.end(); ++bi) {
			Bond *b = *bi;
			connected[b] = (BondOrder) (3 - connected[b]);
		}
		makeAssignments(bonds, connected, curAssign, assignments, allowCharged);
		if (assignments->size() > 0) {
			for (std::vector<Bond *>::iterator bi =
			flip.begin(); bi != flip.end(); ++bi) {
				Bond *b = *bi;
				const Bond::Atoms &bondAtoms =
						b->atoms();
				for (Bond::Atoms::const_iterator
				bai = bondAtoms.begin();
				bai != bondAtoms.end(); ++bai) {
					Atom *a = *bai;
					if (atoms.find(a)
					!= atoms.end())
						continue;
					Element e=a->element();
					if (e == Element::O) {
						if (connected[b] == 1)
							a->setComputedIdatmType("O3");
						else
							a->setComputedIdatmType("O2");
					} else if (e == Element::S) {
						if (connected[b] == 1)
							a->setComputedIdatmType("S3");
						else
							a->setComputedIdatmType("S2");
					}
				}
			}
			break;
		}
		for (std::vector<Bond *>::iterator bi = flip.begin(); bi != flip.end(); 
		++bi) {
			Bond *b = *bi;
			connected[b] = (BondOrder) (3 - connected[b]);
		}
	}
}

}
