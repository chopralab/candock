#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "molecule.hpp"
#include "geom3d/geom3d.hpp"
#include "fragmenter/fragmenter.hpp"
#include "bond.hpp"
#include "idatm.hpp"
using namespace std;

namespace Idatm {
	enum BondOrder { AMBIGUOUS, SINGLE, DOUBLE }; // need SINGLE==1, DOUBLE==2
	void compute_idatm_type(Molib::Molecules &mols) {
		for (auto &molecule : mols)
			compute_idatm_type(molecule);
		dbgmsg(mols);
	}
	void makeAssignments(Molib::BondSet &bonds, map<Molib::Bond*, BondOrder> &connected,
		map<Molib::Bond*, int> &curAssign, vector<map<Molib::Bond*, int>> *assignments, 
		//~ Molib::Bonds&, bool allowCharged=false);
		bool allowCharged=false);
	void flipAssign(Molib::BondVec &flippable, set<Molib::Atom*> &atoms,
		Molib::BondSet &bonds, map<Molib::Bond*, BondOrder> &connected,
		//~ vector<map<Molib::Bond*, int> > *assignments, Molib::Bonds &bond_map, bool allowCharged=false);
		vector<map<Molib::Bond*, int> > *assignments, bool allowCharged=false);
	void uncertainAssign(std::vector<Molib::Atom *> &uncertain,
		std::map<Molib::Atom *, Molib::Bond *> &uncertain2bond,
		std::set<Molib::Bond *> &bonds, std::map<Molib::Bond *, BondOrder> &connected,
		std::vector<std::map<Molib::Bond *, int> > *assignments,
		std::vector<std::vector<Molib::Atom *> > *assignedUncertains,
		//~ Molib::Bonds &bond_map,	bool allowCharged=false);
		bool allowCharged=false);
	//~ bool isN2plus(std::map<Molib::Bond *, int> *bestAssignment, const Molib::Atom &a) {
	//~ bool isN2plus(std::map<Molib::Bond *, int> *bestAssignment, const Molib::BondSet &bonds) {
	bool isN2plus(std::map<Molib::Bond *, int> *bestAssignment, const Molib::BondVec &bonds) {
		int sum = 0, target = 4;
		//~ for (Molib::Atom::Molib::Bonds::const_iterator bi = bonds.begin(); bi != bonds.end();
		//~ ++bi) {
		//~ for (auto &adj_a : a) {
		for (auto &b : bonds) {
			//~ Molib::Bond *b = *bi;
			//~ Molib::Bond *b = bond_map.at((a.atom_number() < adj_a.atom_number() ?
				//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a)));
			std::map<Molib::Bond *, int>::iterator bai = bestAssignment->find(b);
			if (bai == bestAssignment->end())
				target -= 1;
			else
				sum += bai->second;
		}
		return sum == target;
	}

	//~ bool isOarPlus(std::map<Molib::Bond *, int> *bestAssignment, const Molib::Atom &a) {
	//~ bool isOarPlus(std::map<Molib::Bond *, int> *bestAssignment, const Molib::BondSet &bonds) {
	bool isOarPlus(std::map<Molib::Bond *, int> *bestAssignment, const Molib::BondVec &bonds) {
		int sum = 0, target = 3;
		//~ for (Molib::Atom::Molib::Bonds::const_iterator bi = bonds.begin(); bi != bonds.end();
		//~ ++bi) {
		//~ for (auto &adj_a : a) {
		for (auto &b : bonds) {
			//~ Molib::Bond *b = *bi;
			//~ Molib::Bond *b = bond_map.at((a.atom_number() < adj_a.atom_number() ?
				//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a)));
			std::map<Molib::Bond *, int>::iterator bai = bestAssignment->find(b);
			if (bai == bestAssignment->end())
				return false;
			else
				sum += bai->second;
		}
		return sum == target;
	}

	//~ bool isN3plusOkay(const std::vector<Molib::Atom *> &primary) {
	bool isN3plusOkay(const Molib::Atom &a) {
		//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
		//~ bi != primary.end(); ++bi) {
		for (auto &bondee : a) {
		//~ for (auto &bond : a) {
			//~ Molib::Atom &bondee = bond.second_atom();
			//~ Molib::Atom *bondee = *bi;
			//~ Symbol bondeeType = bondee->idatmType();
			//~ const string bondeeType = bondee.idatm_type();
			const string bondeeType = bondee.idatm_type_unmask();
	
			if (bondeeType != "C3" && bondeeType != "H" && bondeeType != "D") {
				return false;
			}
		}
		return true;
	}

	void invertUncertains(std::vector<Molib::Atom *> &uncertain,
		std::map<Molib::Atom *, Molib::Bond *> &uncertain2bond,
		std::map<Molib::Bond *, BondOrder> *connected) {
		for (std::vector<Molib::Atom *>::iterator ui = uncertain.begin();
		ui != uncertain.end(); ++ui) {
			Molib::Atom *a = *ui;
			Molib::Bond *b = uncertain2bond[a];
			(*connected)[b] = (BondOrder) (3 - (*connected)[b]);
			if (a->idatm_type_unmask() == "C3")
				a->set_idatm_type("C2");
			else if (a->idatm_type_unmask() == "C2")
				a->set_idatm_type("C3");
			else if (a->idatm_type_unmask() == "Npl")
				a->set_idatm_type("N2");
			else if (a->idatm_type_unmask() == "N2")
				a->set_idatm_type("Npl");
			else
				std::cerr << "Unknown redo atom type: "
							<< a->idatm_type_unmask() << "\n";
		}
	}

	std::map<Molib::Bond *, int> *findBestAssignment(std::vector<std::map<Molib::Bond *, int> > &assignments,
						//~ std::vector<Molib::Fragmenter::Ring> &systemRings, Molib::Bonds &bond_map) {
						std::vector<Molib::Ring> &systemRings) {
		if (assignments.size() == 0)
			return NULL;
		if (assignments.size() == 1)
			return &assignments[0];
	
		// prefer aromatic if possible (and avoid anti-aromatic)
		std::set<int> okayAssignments;
		std::map<int, int> numPlus;
		int bestAro = 0;
		for (unsigned int i = 0; i < assignments.size(); ++i) {
			std::map<Molib::Bond *, int> &assignment = assignments[i];
			std::map<Molib::Atom *, int> sumOrders, sumBonds;
			for (std::map<Molib::Bond *, int>::iterator ai = assignment.begin();
			ai != assignment.end(); ++ai) {
				Molib::Bond *b = (*ai).first;
				int order = (*ai).second;
				//~ const Molib::Bond::Molib::Atoms &bondAtoms = b->atoms();
				//~ const Molib::Atom::Vec &bondAtoms {&b->first_atom(), &b->second_atom()};
				const Molib::Atom::Vec &bondAtoms {&b->atom1(), &b->atom2()};
				//~ for (Molib::Bond::Molib::Atoms::const_iterator bi =
				//~ bondAtoms.begin(); bi != bondAtoms.end(); ++bi) {
					//~ Molib::Atom *a = *bi;
				for (auto &a : bondAtoms) {
					sumOrders[a] += order;
					sumBonds[a]++;
				}
			}
			// account for N2+
			for (std::map<Molib::Atom *, int>::iterator soi = sumOrders.begin();
			soi != sumOrders.end(); ++soi) {
				Molib::Atom &a = *(*soi).first;
				
				Molib::BondVec primary_bonds;
				//~ for (auto &adj_a : a) {
					//~ primary_bonds.push_back(&bond_map.at((a.atom_number() < adj_a.atom_number() ?
					//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a))));
				//~ }
				//~ for (auto &bond : a) {
					//~ primary_bonds.push_back(&bond);
				//~ }
				for (auto &pbond : a.get_bonds()) {
					primary_bonds.push_back(pbond);
				}
				if (a.element() == Molib::Element::N) {
					//~ if (isN2plus(&assignment, a->primaryBonds()))
					//~ if (isN2plus(&assignment, *a))
					if (isN2plus(&assignment, primary_bonds))
						numPlus[i] += 1;
				} else if (a.element() == Molib::Element::O) {
					//~ if (isOarPlus(&assignment, a->primaryBonds()))
					//~ if (isOarPlus(&assignment, *a))
					if (isOarPlus(&assignment, primary_bonds))
						numPlus[i] += 1;
				}
			}
			int numAro = 0;
			for (std::vector<Molib::Ring>::iterator ri = systemRings.begin();
			ri != systemRings.end(); ++ri) {
				Molib::Ring ring = *ri;
				int piElectrons = 0;
				//~ const Molib::Fragmenter::Ring::Molib::Atoms &atoms = ring.atoms();
				//~ for (Molib::Fragmenter::Ring::Molib::Atoms::const_iterator ai = atoms.begin();
				//~ ai != atoms.end(); ++ai) {
					//~ Molib::Atom *a = *ai;
				for (auto &a : ring) {
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
					//~ } else if (a->primaryBonds().size() == 2 && sum == 2) {
					} else if (a->size() == 2 && sum == 2) {
						piElectrons++;
					} else {
						// aromatic not possible with this assignment
						piElectrons = 0;
						break;
					}
				}
				if (piElectrons % 4 == 2)
					//~ numAro += atoms.size();
					numAro += ring.size();
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
		std::vector<std::map<Molib::Bond *, int> >::iterator bestAssignment, ai;
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
			std::map<Molib::Bond *, int> &assignment = *ai;
			int orderSum = 0;
			for (std::map<Molib::Bond *, int>::iterator i = assignment.begin();
			i != assignment.end(); ++i) {
				//~ val += (*i).first->sqlength() * (*i).second;
				//~ Molib::Atom &batom1 = (*i).first->first_atom(), 
					//~ &batom2 = (*i).first->second_atom();
				Molib::Atom &batom1 = (*i).first->atom1(), 
					&batom2 = (*i).first->atom2();
				val += batom1.crd().distance_sq(batom2.crd()) * (*i).second;
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

	template <class Item>
	void generatePermutations(std::vector<Item *> &items,
					std::vector<std::vector<Item *> > *permutations) {
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

	void uncertainAssign(std::vector<Molib::Atom *> &uncertain,
		std::map<Molib::Atom *, Molib::Bond *> &uncertain2bond,
		std::set<Molib::Bond *> &bonds, std::map<Molib::Bond *, BondOrder> &connected,
		std::vector<std::map<Molib::Bond *, int> > *assignments,
		std::vector<std::vector<Molib::Atom *> > *assignedUncertains,
		//~ Molib::Bonds &bond_map,	bool allowCharged)
		bool allowCharged) {
		std::map<Molib::Bond *, int> curAssign;
		std::vector<std::vector<Molib::Atom *> > permutations;
		generatePermutations<Molib::Atom>(uncertain, &permutations);
		// if we find an assignment involving changing N atoms, have to
		// try all permutations that involve changing no more than N atoms
		unsigned int limit = uncertain.size();
		for (std::vector<std::vector<Molib::Atom *> >::iterator pi =
		permutations.begin(); pi != permutations.end(); ++pi) {
			std::vector<Molib::Atom *> &uncertain = *pi;
			if (uncertain.size() > limit)
				break;
			invertUncertains(uncertain, uncertain2bond, &connected);
			unsigned int numPrevAssigned = assignments->size();
			//~ makeAssignments(bonds, connected, curAssign, assignments, allowCharged);
			//~ makeAssignments(bonds, connected, curAssign, assignments, 
						//~ bond_map, allowCharged);
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

	void flipAssign(Molib::BondVec &flippable, set<Molib::Atom*> &atoms,
		Molib::BondSet &bonds, map<Molib::Bond*, BondOrder> &connected,
		vector<map<Molib::Bond*, int> > *assignments, 
		//~ Molib::Bonds &bond_map, bool allowCharged) {
		bool allowCharged) {
		map<Molib::Bond*, int> curAssign;
		vector<Molib::BondVec > permutations;
		generatePermutations<Molib::Bond>(flippable, &permutations);
		for (vector<Molib::BondVec>::iterator pi =
		permutations.begin(); pi != permutations.end(); ++pi) {
			Molib::BondVec &flip = *pi;
			for (Molib::BondVec::iterator bi =
			flip.begin(); bi != flip.end(); ++bi) {
				Molib::Bond *b = *bi;
				connected[b] = (BondOrder) (3 - connected[b]);
			}
			//~ makeAssignments(bonds, connected, curAssign, assignments, allowCharged);
			//~ makeAssignments(bonds, connected, curAssign, assignments, 
					//~ bond_map, allowCharged);
			makeAssignments(bonds, connected, curAssign, assignments, 
					allowCharged);
			if (assignments->size() > 0) {
				for (Molib::BondVec::iterator bi =
				flip.begin(); bi != flip.end(); ++bi) {
					Molib::Bond *b = *bi;
					//~ const Molib::Bond::Molib::Atoms &bondAtoms =
							//~ b->atoms();
					//~ const Molib::Atom::Vec &bondAtoms {&b->first_atom(), &b->second_atom()};
					const Molib::Atom::Vec &bondAtoms {&b->atom1(), &b->atom2()};
					//~ for (Molib::Bond::Molib::Atoms::const_iterator
					//~ bai = bondAtoms.begin();
					//~ bai != bondAtoms.end(); ++bai) {
					for (auto &a : bondAtoms) {
						//~ Molib::Atom *a = *bai;
						if (atoms.find(a)
						!= atoms.end())
							continue;
						Molib::Element e=a->element();
						if (e == Molib::Element::O) {
							if (connected[b] == 1)
								//~ a->setComputedIdatmType("O3");
								a->set_idatm_type("O3");
							else
								//~ a->setComputedIdatmType("O2");
								a->set_idatm_type("O2");
						} else if (e == Molib::Element::S) {
							if (connected[b] == 1)
								//~ a->setComputedIdatmType("S3");
								a->set_idatm_type("S3");
							else
								//~ a->setComputedIdatmType("S2");
								a->set_idatm_type("S2");
						}
					}
				}
				break;
			}
			for (Molib::BondVec::iterator bi = flip.begin(); bi != flip.end(); 
			++bi) {
				Molib::Bond *b = *bi;
				connected[b] = (BondOrder) (3 - connected[b]);
			}
		}
	}

	void makeAssignments(Molib::BondSet &bonds, map<Molib::Bond*, BondOrder> &connected,
		map<Molib::Bond*, int> &curAssign, vector<map<Molib::Bond*, int>> *assignments, 
		//~ Molib::Bonds &bond_map, bool allowCharged) {
		bool allowCharged) {
		Molib::Bond *assignTarget = *(bonds.begin());
		bonds.erase(bonds.begin());
		bool assign1okay = true, assign2okay = true;
		// see if this assignment completes the bonds of either connected
		// atom and which assignments work
		//~ const Molib::Atom::Set bondAtoms {&assignTarget->first_atom(), &assignTarget->second_atom()};
		const Molib::Atom::Set bondAtoms {&assignTarget->atom1(), &assignTarget->atom2()};
		//~ const Molib::Bond::Molib::Atoms &bondAtoms = assignTarget->atoms(); 
		//~ for (Molib::Bond::Molib::Atoms::const_iterator ai = bondAtoms.begin();
		//~ ai != bondAtoms.end(); ++ai) {
			//~ Molib::Atom *end = *ai;
		for (auto &end : bondAtoms) {
			bool complete = true;
			int sum = 0;
			//~ const Molib::Atom::Molib::Bonds &atomBonds = end->primaryBonds();
			Molib::BondVec atomBonds;
			//~ for (auto &adj_a : *end)
				//~ atomBonds.push_back(&bond_map.at((end->atom_number() < adj_a.atom_number() 
					//~ ? make_pair(end, &adj_a) : make_pair(&adj_a, end))));
			//~ for (auto &bond : *end)
				//~ atomBonds.push_back(&bond);
			for (auto &pbond : end->get_bonds())
				atomBonds.push_back(pbond);
			// implied proton treated the same as ambiguous non-ring bond
			bool hasAmbiguous = atomBonds.size() == 2;
			//~ for (Molib::Atom::Molib::Bonds::const_iterator bi = atomBonds.begin();
			//~ bi != atomBonds.end(); ++bi) {
				//~ Molib::Bond *b = *bi;
			for (auto &b : atomBonds) {
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
				//~ makeAssignments(bonds, connected, curAssign,
							//~ assignments, allowCharged);
				//~ makeAssignments(bonds, connected, curAssign,
							//~ assignments, bond_map, allowCharged);
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
				//~ makeAssignments(bonds, connected, curAssign,
							//~ assignments, allowCharged);
				//~ makeAssignments(bonds, connected, curAssign,
							//~ assignments, bond_map, allowCharged);
				makeAssignments(bonds, connected, curAssign,
							assignments, allowCharged);
			} else {
				assignments->push_back(curAssign);
			}
			curAssign.erase(assignTarget);
		}
		bonds.insert(assignTarget);
	}

	bool aromaticGeometry(const Molib::Atom::Set &r) {
		// algorithm from:
		//	Crystallographic Studies of Inter- and Intramolecular 
		//	   Interactions Reflected in Aromatic Character of pi-Electron
		//	   Systems
		//	J. Chem. Inf. Comput. Sci., 1993, 33, 70-78
		double sum = 0.0;
		int bonds = 0;
		//~ set<Molib::BondKey> ring_bonds;		
		//~ for (auto &a : r) {
			//~ for (auto &adj_a : *a) {
				//~ if (r.count(&adj_a) && !ring_bonds.count(Molib::BondKey(&adj_a, a))) {
					//~ ring_bonds.insert(Molib::BondKey(a, &adj_a));
					//~ dbgmsg("aromaticity test : ring bond between atom " << *a << " and atom " << adj_a);
				//~ }
			//~ }
		//~ }
		auto ring_bonds = get_bonds_in(r);

		//~ for (Molib::Fragmenter::Ring::Molib::Bonds::const_iterator bi = r.bonds().begin();
		//~ bi != r.bonds().end(); ++bi) {
			//~ Molib::Bond *b = *bi;
		for (auto &pbond : ring_bonds) {
			Molib::Bond &bond = *pbond;
			const Molib::Element &e1 = bond.atom1().element();
			const Molib::Element &e2 = bond.atom2().element();
			const Geom3D::Coordinate &c1 = bond.atom1().crd();
			const Geom3D::Coordinate &c2 = bond.atom2().crd();
			double d = c1.distance(c2), delta;
			//~ const Molib::Element &e1 = b->atoms()[0]->element();
			//~ const Molib::Element &e2 = b->atoms()[1]->element();
			//~ Coord c1 = b->atoms()[0]->coord();
			//~ Coord c2 = b->atoms()[1]->coord();
			//~ double d = distance(c1, c2), delta;
	
			if (e1 == Molib::Element::C && e2 == Molib::Element::C) {
				delta = d - 1.38586;
			} else if ((e1 == Molib::Element::C || e2 == Molib::Element::C) &&
			  (e1 == Molib::Element::N || e2 == Molib::Element::N)) {
				delta = d - 1.34148;
			} else
				continue;
			bonds++;
			sum += delta * delta;
	
		}
		if (bonds == 0)
			return false;
		
		double homa = 1.0 - (792.0 * sum / bonds);
	
		if (homa >= 0.5) {
			dbgmsg("ring is aromatic");
			// aromatic
			return true;
		} else if (bonds * homa < -35.0)
			return false;
	
		return true;
	}

	int freeOxygens(const Molib::Atom &a, map<const Molib::Atom*, int> &heavys, 
		const bool noHyds=false) {
		int freeOxygens = 0;
		//~ for (auto &bond : a) {
		for (auto &bondee : a) {
			//~ Molib::Atom &bondee = bond.second_atom();
			if (bondee.element() == Molib::Element::O) {
				if (heavys[&bondee] != 1)
					continue;
				if (noHyds && bondee.size() != 1)
					continue;
				freeOxygens++;
			}
		}
		return freeOxygens;
	}
//~ int
//~ freeOxygens(std::vector<Molib::Atom *> &primary, std::map<Molib::Atom *, int> &heavys, bool noHyds)
//~ {
	//~ int freeOxygens = 0;
	//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
	  //~ bi != primary.end(); ++bi) {
		//~ Molib::Atom *bondee = *bi;
		//~ if (bondee->element() == Molib::Element::O) {
			//~ if (heavys[bondee] != 1)
				//~ continue;
			//~ if (noHyds && bondee->primaryNeighbors().size() != 1)
				//~ continue;
			//~ freeOxygens++;
		//~ }
	//~ }
//~ 
	//~ return freeOxygens;
//~ }
	
	
	void compute_idatm_type(Molib::Molecule &molecule) {
		// angle values used to discriminate between hybridization states
		const double angle23val1 = 115.0;
		const double angle23val1_tmax = 116.5;
		const double angle23val1_tmin = 113.5;
		const double angle23val2 = 122.0;
		const double angle12val = 160.0;
	
		// bond length cutoffs from hybridization discrimination
		// p3... = pass 3 cutoffs; p4... = pass 4 cutoffs
		const double p3c1c1 = 1.22 * 1.22;
		const double p3c2c = 1.41 * 1.41;
		const double p3c2n = 1.37 * 1.37;
		const double p3n1c1 = 1.20 * 1.20;
		const double p3n3c = 1.38 * 1.38;
		const double p3n3n3 = 1.43 * 1.43;
		const double p3n3n2 = 1.41 * 1.41;
		const double p3n1o1 = 1.21 * 1.21;
		const double p3c1o1 = 1.17 * 1.17;
		const double p3o2c2 = 1.30 * 1.30;
		const double p3o2as = 1.685 * 1.685;
		const double p3o2o3 = 1.338 * 1.338;
		const double p3s2c2 = 1.76 * 1.76;
		const double p3s2as = 2.11 * 2.11;
		const double p4c3c = 1.53 * 1.53;
		const double p4c3n = 1.46 * 1.46;
		const double p4c3o = 1.44 * 1.44;
		const double p4n2c = 1.38 * 1.38;
		const double p4n2n = 1.32 * 1.32;
		const double p4c2c = 1.42 * 1.42;
		const double p4c2n = 1.41 * 1.41;
		const double p4ccnd = 1.45 * 1.45;
	
		const double p7cn2nh = 1.3629;
		const double p7nn2nh = 1.3337;
		const double p7on2nh = 1.3485;
	
		// algorithm based on E.C. Meng / R.A. Lewis paper 
		// "Determination of Molecular Topology and Molib::Atomic Hybridization
		// States from Heavy Molib::Atom Coordinates", J. Comp. Chem., v12#7, 891-898
		// and on example code from idatm.f implementation by E.C. Meng
	
		// differences: No boron types.  Double-bonded Npls are split off
		//   as N2.  Sox split into Sxd (sulfoxide), and Son (sulfone).
		//   Carbons in aromatic rings are type Car.  Aromatic oxygens are Oar/Oar+.
		//   Negatively charged oxygens are O2- (planar) and O3- (tetrahedral)
		//   instead of O-.  Sp nitrogens bonded to two atoms are N1+.
		
		//~ const Molib::Atom::IdatmInfoMap &infoMap = Molib::Atom::getIdatmInfoMap();
	
		size_t numAtoms = 0;
		// initialize idatm type in Molib::Atoms
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
			a.set_idatm_type(a.element().name());
			++numAtoms;
		}
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
			//~ a->setComputedIdatmType(a->element().name());
		//~ }
	
		//~ // if molecule is diamond/nanotube, skip atom typing since the
		//~ // ring finding will take forever
		//~ size_t numBonds = mol->numBonds();
		//~ size_t numAtoms = mol->numAtoms();
		//~ if (numBonds - numAtoms > 100 && numBonds / (float) numAtoms > 1.25)
			//~ return;
	
		map<const Molib::Atom*, int> heavys; // number of heavy atoms bonded
		//~ std::map<Molib::Atom *, int> heavys; // number of heavy atoms bonded
		//~ std::vector<Molib::Atom *> primary; // since primaryNeighbors() returns a copy
		size_t hassigned = 0;
	
		// "pass 1":  type hydrogens / deuteriums and compute number of
		// heavy atoms connected to each atom
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
			const Molib::Element &element = a.element();
			dbgmsg(a.element());
			//~ const Molib::Element &element = a->element();
			//~ primary = a->primaryNeighbors();
	
			if (element.number() == 1) {
				// sort out if it's a hydrogen or deuterium
				bool isHyd = true;
				if (a.atom_name()[0] == 'D' || a.atom_name()[0] == 'd')
					isHyd = false;
				//~ for (const char *c = a->name().c_str(); *c != '\0';
				  //~ ++c) {
					//~ if (isalpha(*c)) {
						//~ if (*c == 'd' || *c == 'D') {
							//~ isHyd = false;
						//~ }
						//~ break;
					//~ }
				//~ }
				bool bondedToCarbon = false;
				for (auto &bondee : a) {
				//~ for (auto &bond : a) {
					//~ auto &bondee = bond.second_atom();
					if (bondee.element() == Molib::Element::C) {
						bondedToCarbon = true;
						break;
					}
				}
				//~ bool bondedToCarbon = false;
				//~ for (std::vector<Molib::Atom *>::const_iterator bi =
				  //~ primary.begin(); bi != primary.end(); ++bi) {
				  	//~ Molib::Atom *bondee = *bi;
					//~ if (bondee->element() == Molib::Element::C) {
						//~ bondedToCarbon = true;
						//~ break;
					//~ }
				//~ }
				a.set_idatm_type(bondedToCarbon ? (isHyd ? "HC" : "DC") : (isHyd ? "H" : "D"));
				dbgmsg("hydrogens" << boolalpha << (bondedToCarbon ? (isHyd ? "HC" : "DC") : (isHyd ? "H" : "D")));
				//~ a->setComputedIdatmType(bondedToCarbon ? (isHyd ?
				  //~ "HC" : "DC") : (isHyd ? "H" : "D"));
				hassigned += 1;
			}
	
			int heavyCount = 0;
			for (auto &bondee : a) {
			//~ for (auto &bond : a) {
				//~ auto &bondee = bond.second_atom();
				if (bondee.element().number() > 1) {
					heavyCount++;
				}
			}
			//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
			  //~ bi != primary.end(); ++bi) {
				//~ Molib::Atom *bondee = *bi;
				//~ if (bondee->element().number() > 1) {
					//~ heavyCount++;
				//~ }
			//~ }
			heavys[&a] = heavyCount;
			dbgmsg("pass 1  : " << a.idatm_type() << " heavys[" << a.atom_name() << "] = " << heavys[&a]);
		}
	
		// "pass 1.5": use templates for "infallible" typing of standard
		// residue types
		map<const Molib::Atom*, bool> mapped;
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model) {
			//~ Molib::Atom &first_atom = chain.first().first(); // take care of the first and the last atom
			//~ if (first_atom.atom_name() == "N") { first_atom.set_idatm_type("N3+"); mapped[&first_atom] = true; }
			for (auto &residue : chain) {
				auto it = help::standard_residues.find(residue.resn());
				if (it != help::standard_residues.end())
					for (auto &a : residue) {
						if (!mapped[&a]) {
							auto it2 = it->second.find(a.atom_name());
							if (it2 == it->second.end()) throw Error("die : cannot find atom name " + a.atom_name() + " of template residue");
							a.set_idatm_type(it2->second);
							mapped[&a] = true;
						}
						dbgmsg("pass 1.5  : " << a.idatm_type() << " " << a.atom_name() << " " << a.atom_number());
				}
			}
		}
		//~ std::map<const Molib::Atom *, bool> mapped;
		//~ const Molecule::Residues &res = mol->residues();
		//~ for (Molecule::Residues::const_iterator ri = res.begin();
		     //~ ri != res.end(); ++ri) {
			//~ Residue *r = *ri;
			//~ try {
				//~ std::vector<Molib::Atom *> templated = r->templateAssign(
				 //~ &Molib::Atom::setComputedIdatmType, "idatm", "templates", "idatmres");
				//~ for (std::vector<Molib::Atom *>::iterator ai =
				//~ templated.begin(); ai != templated.end(); ++ai)
					//~ mapped[*ai] = true;
			//~ } catch (TA_NoTemplate) {
				//~ // don't care
			//~ } catch (...) {
				//~ throw;
			//~ }
		//~ }
	
		/*
		std::cerr << "  computeIdatmTypes() hydrogens " << hassigned;
		std::cerr << ", template atoms " << mapped.size();
		std::cerr << ", unassigned " << numAtoms - (hassigned + mapped.size()) << std::endl;
		*/
	
		if (hassigned + mapped.size() == numAtoms)
			return;		// All atoms assigned.
	
		// "pass 2": elements that are typed only by element type
		// and valences > 1
		map<Molib::Atom*, int> redo;
		//~ std::map<Molib::Atom *, int> redo;
		//~ std::set<Molib::Atom *> ambiguousVal2Cs;
		set<Molib::Atom*> ambiguousVal2Cs;
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
			if (mapped[&a])
				continue;
			const Molib::Element &element = a.element();
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
			//~ if (mapped[a])
				//~ continue;
			//~ const Molib::Element &element = a->element();
	
			// undifferentiated types
			if (element >= Molib::Element::He && element <= Molib::Element::Be
			  || element >= Molib::Element::Ne && element <= Molib::Element::Si
			  || element >= Molib::Element::Cl) {
			  	a.set_idatm_type(element.name());
			  	dbgmsg(a.idatm_type());
			  	//~ a->setComputedIdatmType(element.name());
				continue;
			}
	
			// valence 4
			//	C must be sp3 (C3)
			//	N must be part of a quaternary amine (N3+) 
			//	P must be part of a phosphate (Pac), a P-oxide (Pox)
			//		or a quaternary phosphine (P3+)
			//	S must be part of a sulfate, sulfonate or sulfamate
			//		(Sac), or sulfone (Son)
			//~ primary = a->primaryNeighbors();
			if (a.size() == 4) {
			//~ if (primary.size() == 4) {
				if (element == Molib::Element::C) {
					a.set_idatm_type("C3");
					//~ a->setComputedIdatmType("C3");
				} else if (element == Molib::Element::N) {
					a.set_idatm_type("N3+");
					//~ a->setComputedIdatmType("N3+");
				} else if (element == Molib::Element::P) {
					int freeOxys = freeOxygens(a, heavys);
					//~ int freeOxys = freeOxygens(primary, heavys);
					if (freeOxys >= 2)
						a.set_idatm_type("Pac");
						//~ a->setComputedIdatmType("Pac");
					else if (freeOxys == 1)
						a.set_idatm_type("Pox");
						//~ a->setComputedIdatmType("Pox");
					else
						a.set_idatm_type("P3+");
						//~ a->setComputedIdatmType("P3+");
				} else if (element == Molib::Element::S) {
					int freeOxys = freeOxygens(a, heavys);
					//~ int freeOxys = freeOxygens(primary, heavys);
					if (freeOxys >= 3) {
						a.set_idatm_type("Sac");
					} else if (freeOxys >= 1) {
						a.set_idatm_type("Son");
					} else {
						a.set_idatm_type("S");
					}
					//~ if (freeOxys >= 3) {
						//~ a->setComputedIdatmType("Sac");
					//~ } else if (freeOxys >= 1) {
						//~ a->setComputedIdatmType("Son");
					//~ } else {
						//~ a->setComputedIdatmType("S");
					//~ }
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
			else if (a.size() == 3) {
				double avgAngle = 0.0;
				for (int n1 = 0; n1 < 3; ++n1) {
					for (int n2 = n1 + 1; n2 < 3; ++n2) {
						avgAngle += Geom3D::degrees(Geom3D::angle(
							a[n1].crd(),
							a.crd(),
							a[n2].crd()));
					}
				}
				//~ for (auto it1 = a.begin(); it1 != a.end(); ++it1) {
					//~ auto it2 = it1;
					//~ for (++it2; it2 != a.end(); ++it2) {
						//~ auto &a1 = it1->second_atom();
						//~ auto &a2 = it2->second_atom();
						//~ avgAngle += Geom3D::degrees(Geom3D::angle(
							//~ a1.crd(),
							//~ a.crd(),
							//~ a2.crd()));
					//~ }
				//~ }
				avgAngle /= 3.0;
			//~ else if (primary.size() == 3) {
				//~ Real avgAngle = 0.0;
				//~ for (int n1 = 0; n1 < 3; ++n1) {
					//~ for (int n2 = n1 + 1; n2 < 3; ++n2) {
						//~ avgAngle += angle(
						  //~ primary[n1]->coord(), a->coord(),
						  //~ primary[n2]->coord());
					//~ }
				//~ }
				//~ avgAngle /= 3.0;
	
				if (element == Molib::Element::C) {
					bool c3 = false;
					if (avgAngle < angle23val1_tmin)
						c3 = true;
					else if (avgAngle < angle23val1_tmax) {
						double minSqDist = -1.0;
						//~ Real minSqDist = -1.0;
						for (int n1 = 0; n1 < 3; ++n1) {
							double sqd = a.crd().distance_sq(a[n1].crd());
						//~ for (auto &bond : a) {
							//~ auto &bondee = bond.second_atom();
							//~ double sqd = a.crd().distance_sq(bondee.crd());
							//~ Real sqd = a->coord().sqdistance(primary[n1]->coord());
							if (minSqDist < 0.0 || sqd < minSqDist)
								minSqDist = sqd;
						}
						if (minSqDist > p4c3c)
							c3 = true;
						else if (minSqDist > p4c2c && avgAngle < angle23val1)
							c3 = true;
					}
					if (c3)
						a.set_idatm_type("C3");
						//~ a->setComputedIdatmType("C3");
					else
						a.set_idatm_type(freeOxygens(a, heavys, true) >= 2 ?  "Cac" : "C2");
						//~ a->setComputedIdatmType(freeOxygens(
							//~ primary, heavys, true) >= 2
							//~ ? "Cac" : "C2");
				} else if (element == Molib::Element::N) {
					if (avgAngle < angle23val1)
						a.set_idatm_type("N3");
						//~ a->setComputedIdatmType("N3");
					else
						a.set_idatm_type(freeOxygens(a, heavys) >= 2 ?  "Ntr" : "Npl");
						//~ a->setComputedIdatmType(freeOxygens(
						  //~ primary, heavys) >= 2 ? "Ntr":"Npl");
				} else if (element == Molib::Element::S) {
					bool hasOxy = false;
					for (int i = 0; i < 3; ++i) {
						if (a[i].element() == Molib::Element::O) {
					//~ for (auto &bond : a) {
						//~ auto &bondee = bond.second_atom();
						//~ if (bondee.element() == Molib::Element::O) {
						//~ if (primary[i]->element() ==
						  //~ Molib::Element::O) {
						  	hasOxy = true;
							break;
						}
					}
					a.set_idatm_type(hasOxy ? "Sxd" : "S3+");
					//~ a->setComputedIdatmType(hasOxy ? "Sxd" : "S3+");
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
			else if (a.size() == 2) {
				double ang = Geom3D::degrees(Geom3D::angle(a[0].crd(), a.crd(), a[1].crd()));
				//~ double ang = Geom3D::degrees(Geom3D::angle(a.first().second_atom().crd(), 
					//~ a.crd(), a.last().second_atom().crd()));
			//~ else if (primary.size() == 2) {
				//~ Point coord[2];
				//~ int coordInd = 0;
				//~ for (std::vector<Molib::Atom *>::const_iterator bi =
				  //~ primary.begin(); bi != primary.end(); ++bi) {
				  	//~ Molib::Atom *other = *bi;
					//~ coord[coordInd++] = other->coord();
				//~ }
				//~ Real ang = angle(coord[0], a->coord(), coord[1]);
	
				if (element == Molib::Element::C) {
					if (ang < angle23val1) {
						a.set_idatm_type("C3");
						redo[&a] = 1;
						//~ a->setComputedIdatmType("C3");
						//~ redo[a] = 1;
						if (ang > angle23val1_tmin)
							ambiguousVal2Cs.insert(&a);
					} else if (ang < angle12val) {
						a.set_idatm_type("C2");
						//~ a->setComputedIdatmType("C2");
						if (ang < angle23val2) {
							redo[&a] = 3;
						} else {
							// allow ring bond-order code
							// to change this assignment
							redo[&a] = -1;
						}
						if (ang < angle23val1_tmax)
							ambiguousVal2Cs.insert(&a);
					} else {
						a.set_idatm_type("C1");
						//~ a->setComputedIdatmType("C1");
					}
				} else if (element == Molib::Element::N) {
					if (ang < angle23val1) {
						a.set_idatm_type("N3");
						//~ a->setComputedIdatmType("N3");
						redo[&a] = 2;
					} else {
						a.set_idatm_type(ang < angle12val ? "Npl" : "N1+");
						//~ a->setComputedIdatmType(
						  //~ ang < angle12val ?  "Npl" : "N1+");
					}
				} else if (element == Molib::Element::O) {
					a.set_idatm_type("O3");
					//~ a->setComputedIdatmType("O3");
				} else if (element == Molib::Element::S) {
					a.set_idatm_type("S3");
					//~ a->setComputedIdatmType("S3");
				}
			}
			dbgmsg("pass 2  : " << a.idatm_type());
		}
	
		// "pass 3": determine types of valence 1 atoms.  These were typed
		// by element only in previous pass, but can be typed more accurately
		// now that the atoms they are bonded to have been typed.  Molib::Bond
		// lengths are used in this pass.  
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
			if (a.size() != 1)
				continue;
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
	//~ 
			//~ primary = a->primaryNeighbors();
			//~ if (primary.size() != 1)
				//~ continue;
			
			Molib::Atom &bondee = a.first();
			//~ Molib::Atom &bondee = a.first().second_atom();
			double sqlen = bondee.crd().distance_sq(a.crd());
			const string bondeeType = bondee.idatm_type_unmask();
			//~ Molib::Atom *bondee = *(primary.begin());
			//~ Real sqlen = sqdistance(bondee->coord(), a->coord());
			//~ Symbol bondeeType = bondee->idatmType();
	
			
			if (a.idatm_type_unmask() == "C") {
				if (mapped[&a])
					continue;
			//~ if (a->idatmType() == "C") {
				//~ if (mapped[a])
					//~ continue;
				if (sqlen <= p3c1c1 && bondeeType == "C1") {
					a.set_idatm_type("C1");
					//~ a->setComputedIdatmType("C1");
				} else if (sqlen <= p3c2c &&
				  bondee.element() == Molib::Element::C) {
					a.set_idatm_type("C2");
					//~ a->setComputedIdatmType("C2");
				} else if (sqlen <= p3c2n &&
				  bondee.element() == Molib::Element::N) {
					a.set_idatm_type("C2");
					//~ a->setComputedIdatmType("C2");
				} else if (sqlen <= p3c1o1 &&
				  bondee.element() == Molib::Element::O &&
				  //~ bondee->primaryNeighbors().size() == 1) {
				  bondee.size() == 1) {
				  	a.set_idatm_type("C1-");
				  	//~ a->setComputedIdatmType("C1-");
				} else {
					a.set_idatm_type("C3");
					//~ a->setComputedIdatmType("C3");
				}
			} else if (a.idatm_type_unmask() == "N") {
			//~ } else if (a->idatmType() == "N") {
				if (mapped[&a])
					continue;
				//~ if (mapped[a])
					//~ continue;
				if ((sqlen <= p3n1c1 && bondeeType == "C1" ||
				  bondeeType == "N1+") || (sqlen < p3n1o1 &&
				  //~ bondee->element() == Molib::Element::O)) {
				  bondee.element() == Molib::Element::O)) {
					a.set_idatm_type("N1");
					//~ a->setComputedIdatmType("N1");
				} else if (sqlen > p3n3c &&
				  (bondeeType == "C2" || bondeeType == "C3")) {
					a.set_idatm_type("N3");
				//~ } else if (sqlen > p3n3c &&
				  //~ (bondeeType == "C2" || bondeeType == "C3")) {
					//~ a->setComputedIdatmType("N3");
				} else if ((sqlen > p3n3n3 && bondeeType == "N3") || 
					(sqlen > p3n3n2 && bondeeType == "Npl")) {
					a.set_idatm_type("N3");
				//~ } else if ((sqlen > p3n3n3 && bondeeType == "N3") ||
				  //~ (sqlen > p3n3n2 && bondeeType == "Npl")) {
					//~ a->setComputedIdatmType("N3");
				} else if (bondee.element() == Molib::Element::C ||
				  bondee.element() == Molib::Element::N) {
					a.set_idatm_type("Npl");
					//~ a->setComputedIdatmType("Npl");
				} else {
					a.set_idatm_type("N3");
					//~ a->setComputedIdatmType("N3");
				}
			} else if (a.idatm_type_unmask() == "O") {
			//~ } else if (a->idatmType() == "O") {
				if (bondeeType == "Cac" || bondeeType == "Ntr" ||
				  bondeeType == "N1+") {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("O2-");
						a.set_idatm_type("O2-");
				} else if (bondeeType == "Pac" || bondeeType == "Sac"
				  || bondeeType == "N3+" || bondeeType == "Pox"
				  || bondeeType == "Son" || bondeeType == "Sxd") {
					if (mapped[&a])
						continue;
					a.set_idatm_type("O3-");
					//~ a->setComputedIdatmType("O3-");
	
					// pKa of 3rd phosphate oxygen is 7...
					if (bondeeType != "Pac")
						continue;
					//~ std::vector<Molib::Atom *> bondeePrimary =
							//~ bondee->primaryNeighbors();
					int oxys = 0;
					for (auto &bp : bondee) {
					//~ for (auto &bond : bondee) {
						//~ auto &bp = bond.second_atom();
					//~ for (std::vector<Molib::Atom *>::const_iterator bpi =
					//~ bondeePrimary.begin(); bpi !=
					//~ bondeePrimary.end(); ++bpi) {
						//~ Molib::Atom *bp = *bpi;
						if (bp.element() == Molib::Element::O
						&& bp.size() == 1)
						//~ if (bp->element() == Molib::Element::O
						//~ && bp->primaryNeighbors().size() == 1)
							oxys += 1;
					}
					if (oxys < 3)
						continue;
					// if this bond is 0.05 A longer than
					// the other P-O bonds, assume OH
					//~ Real len = distance(bondee->coord(),
									//~ a->coord());
					double len = bondee.crd().distance(a.crd());
					bool longer = true;
					for (auto &bp : bondee) {
					//~ for (auto &bond : bondee) {
						//~ auto &bp = bond.second_atom();
					//~ for (std::vector<Molib::Atom *>::const_iterator bpi =
					//~ bondeePrimary.begin(); bpi !=
					//~ bondeePrimary.end(); ++bpi) {
						//~ Molib::Atom *bp = *bpi;
						if (&bp == &a
						|| bp.size() > 1
						//~ || bp->primaryNeighbors().size() > 1
						|| bp.element() != Molib::Element::O)
							continue;
						//~ if (len < distance(bondee->coord(),
						//~ bp->coord()) + 0.05) {
						if (len < bondee.crd().distance(bp.crd()) + 0.05) {
							longer = false;
							break;
						}
					}
					if (longer)
						//~ a->setComputedIdatmType("O3");
						a.set_idatm_type("O3");
				} else if (sqlen <= p3c1o1 &&
				  bondee.element() == Molib::Element::C &&
				  //~ bondee->primaryNeighbors().size() == 1) {
				  bondee.size() == 1) {
				  	if (!mapped[&a])
						//~ a->setComputedIdatmType("O1+");
						a.set_idatm_type("O1+");
				} else if (sqlen <= p3o2c2 &&
				  bondee.element() == Molib::Element::C) {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("O2");
						a.set_idatm_type("O2");
					if (!mapped[&bondee])
						//~ bondee->setComputedIdatmType("C2");
						bondee.set_idatm_type("C2");
					redo[&bondee] = 0;
				} else if (sqlen <= p3o2as &&
				  bondee.element() == Molib::Element::As) {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("O2");
						a.set_idatm_type("O2");
				} else if (sqlen <= p3o2o3 &&
				  bondee.element() == Molib::Element::O &&
				  //~ bondee->primaryNeighbors().size() == 1) {
				  bondee.size() == 1) {
					// distinguish oxygen molecule from
					// hydrogen peroxide
					if (!mapped[&a])
						//~ a->setComputedIdatmType("O2");
						a.set_idatm_type("O2");
				} else if (sqlen <= p3n1o1 &&
				  bondee.element() == Molib::Element::N &&
				  //~ bondee->primaryNeighbors().size() == 1) {
				  bondee.size() == 1) {
				  	if (!mapped[&a])
						//~ a->setComputedIdatmType("O1");
						a.set_idatm_type("O1");
				} else {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("O3");
						a.set_idatm_type("O3");
				}
			//~ } else if (a->idatmType() == "S") {
			} else if (a.idatm_type_unmask() == "S") {
				if (bondee.element() == Molib::Element::P) {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("S3-");
						a.set_idatm_type("S3-");
				} else if (bondeeType == "N1+") {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("S2");
						a.set_idatm_type("S2");
				} else if (sqlen <= p3s2c2 &&
				  bondee.element() == Molib::Element::C) {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("S2");
						a.set_idatm_type("S2");
					if (!mapped[&bondee])
						//~ bondee->setComputedIdatmType("C2");
						bondee.set_idatm_type("C2");
					redo[&bondee] = 0;
				} else if (sqlen <= p3s2as &&
				  //~ bondee->element() == Molib::Element::As) {
				  bondee.element() == Molib::Element::As) {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("S2");
						a.set_idatm_type("S2");
				} else {
					if (!mapped[&a])
						//~ a->setComputedIdatmType("S3");
						a.set_idatm_type("S3");
				}
			}
			dbgmsg("pass 3  : " << a.idatm_type_unmask());
		}
	
		// "pass 4": re-examine all atoms with non-zero 'redo' values and
		//   retype them if necessary
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
			if (mapped[&a])
				redo[&a] = 0;
	
			if (redo[&a] == 0)
				continue;
			
			bool c3able = false;
			//~ primary = a->primaryNeighbors();
			for (auto &bondee : a) {
			//~ for (auto &bond : a) {
				//~ auto &bondee = bond.second_atom();
			//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
			  //~ bi != primary.end(); ++bi) {
				//~ Molib::Atom *bondee = *bi;
				const Molib::Element &bondeeElement = bondee.element();
				//~ Real sqlen = sqdistance(bondee->coord(),
									//~ a->coord());
				double sqlen = bondee.crd().distance_sq(a.crd());
	
				if (redo[&a] == 1) {
					if ((sqlen <= p4c2c && bondeeElement ==
					  Molib::Element::C) || (sqlen <= p4c2n &&
					  bondeeElement == Molib::Element::N)) {
						a.set_idatm_type("C2");
						//~ a->setComputedIdatmType("C2");
						break;
					}
					if ((sqlen > p4c3c && bondeeElement ==
					  Molib::Element::C) || (sqlen > p4c3n &&
					  bondeeElement == Molib::Element::N) || (sqlen >
					  p4c3o && bondeeElement == Molib::Element::O)) {
						//~ a->setComputedIdatmType("C3");
						a.set_idatm_type("C3");
					}
				} else if (redo[&a] == 2) {
					if ((sqlen <= p4n2c && bondeeElement ==
					  Molib::Element::C) || (sqlen <= p4n2n &&
					  bondeeElement == Molib::Element::N)) {
						// explicit hydrogen(s): N2
						if (heavys[&a] < 2)
							//~ a->setComputedIdatmType("N2");
							a.set_idatm_type("N2");
						else
							//~ a->setComputedIdatmType("Npl");
							a.set_idatm_type("Npl");
						break;
					}
				} else {
					if ((sqlen <= p4c2c && bondeeElement ==
					  Molib::Element::C) || (sqlen <= p4c2n &&
					  bondeeElement == Molib::Element::N)) {
						//~ a->setComputedIdatmType("C2");
						a.set_idatm_type("C2");
						c3able = false;
						break;
					}
					if ((sqlen > p4c3c && bondeeElement ==
					  Molib::Element::C) || (sqlen > p4c3n &&
					  bondeeElement == Molib::Element::N) || (sqlen >
					  p4c3o && bondeeElement == Molib::Element::O)) {
						c3able = true;
					}
					if (sqlen > p4ccnd &&
					  bondeeElement == Molib::Element::C) {
						c3able = true;
					}
				}
			}
			if (c3able)
				//~ a->setComputedIdatmType("C3");
				a.set_idatm_type("C3");
			dbgmsg("pass 4  : " << a.idatm_type_unmask());
		}
	
		// "pass 4.5":  this pass is not in the IDATM paper but is a suggested
		//    improvement mentioned on page 897 of the paper:  find aromatic
		//    ring types.  The method is to:
		//
		//	1) Find all intraresidue rings (actually computed before pass 6)
		//	2) Check that all the atoms of the ring are planar types
		//	3) Check bond lengths around the ring; see if they are
		//		consistent with aromatic bond lengths
		set<Molib::Atom*> ringAssignedNs;
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		if (!help::standard_residues.count(residue.resn())) {
			//~ Molib::Bonds bond_map = Molib::create_bonds(Molib::create_graph(residue)); // ?????????????????????????
			for (auto &a : residue) {
				//~ std::set<const Residue *> mappedResidues;
				//~ for (std::map<const Molib::Atom *, bool>::const_iterator mi = mapped.begin();
						//~ mi != mapped.end(); ++mi) {
					//~ if ((*mi).second)
						//~ mappedResidues.insert((*mi).first->residue());
				//~ }
				//~ std::set<const Molib::Atom *> considerMapped;
				//~ for (std::set<const Residue *>::const_iterator mri = mappedResidues.begin();
						//~ mri != mappedResidues.end(); ++mri) {
					//~ const Residue *r = *mri;
					//~ const Residue::Molib::Atoms &rMolib::Atoms = r->atoms();
					//~ for (Residue::Molib::Atoms::const_iterator rai = rMolib::Atoms.begin(); rai != rMolib::Atoms.end(); ++rai) {
						//~ considerMapped.insert(*rai);
					//~ }
				//~ }
				//~ int ringLimit = 3;
				//~ Molib::Fragmenter::Rings rs = Molib::Fragmenter(residue).identify_rings();
				Molib::Rings rs = Molib::Fragmenter(residue.get_atoms()).identify_rings();
				//~ Molib::Fragmenter::Rings rs = mol->rings(false, ringLimit, &considerMapped);
				// std::cerr << "  computeIdatmTypes() size 3 rings " << rs.size() << std::endl;
				//~ if (rs.size() < 20) {
					//~ // not something crazy like an averaged structure...
					//~ ringLimit = 6;
					//~ rs = mol->rings(false, ringLimit, &considerMapped);
					//~ // std::cerr << "  computeIdatmTypes() size 6 rings " << rs.size() << std::endl;
					//~ if (rs.size() < 20) {
						//~ // not something crazy like a nanotube...
						//~ ringLimit = 0;
						//~ rs = mol->rings(false, ringLimit, &considerMapped);
					//~ }
				//~ }
				// screen out rings with definite non-planar types
				Molib::Rings planarRings;
				for (auto &r : rs) {
				//~ for (Molib::Fragmenter::Rings::const_iterator ri = rs.begin(); ri != rs.end(); ++ri) {
					//~ const Molib::Fragmenter::Ring &r = *ri;
			
					//~ if (r.atoms().size() == 3) {
					if (r.size() == 3) {
						for (auto &pa : r) {
							Molib::Atom &a = *pa;
						//~ for (Molib::Fragmenter::Ring::Molib::Atoms::const_iterator ai = r.atoms().begin();
						//~ ai != r.atoms().end(); ++ai) {
							//~ Molib::Atom *a = *ai;
							//~ if (a->element() == Molib::Element::C)
							if (a.element() == Molib::Element::C)
								//~ a->setComputedIdatmType("C3");
								a.set_idatm_type("C3");
						}
						continue;
					}
					bool planarTypes = true;
					bool allMapped = true;
					bool allPlanar = true;
					int numOxygens = 0;
					for (auto &pa : r) {
						Molib::Atom &a = *pa;
					//~ for (Molib::Fragmenter::Ring::Molib::Atoms::const_iterator ai = r.atoms().begin();
					//~ ai != r.atoms().end(); ++ai) {
						//~ const Molib::Atom *a = *ai;
						//~ Symbol idatmType = a->idatmType();
						const string idatmType = a.idatm_type_unmask();
						if (a.element() == Molib::Element::O)
							numOxygens++;
			
						//~ primary = a->primaryNeighbors();
						//~ if (primary.size() > 3) {
						if (a.size() > 3) {
							allPlanar = planarTypes = false;
							break;
						}
			
						if (idatmType != "C2" && idatmType != "Npl" &&
						  idatmType != "Sar" && idatmType != "O3" &&
						  idatmType != "S3" && idatmType != "N3" &&
						  idatmType != "Oar" && idatmType != "Oar+" && idatmType != "P" &&
						  idatmType != "Car" && idatmType != "N2" && idatmType != "N2+" &&
						  //~ !(idatmType == "C3" && primary.size()==2)) {
						  !(idatmType == "C3" && a.size()==2)) {
							allPlanar = planarTypes = false;
							break;
						} else if (idatmType == "O3" || idatmType == "S3" ||
						  idatmType == "N3" || idatmType == "C3") {
						  	allPlanar = false;
						}
			
						//~ if (mapped[a]) {
						if (mapped[&a]) {
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
					//~ if (r.atoms().size() == 5 && numOxygens > 1 && numOxygens < 5)
					if (r.size() == 5 && numOxygens > 1 && numOxygens < 5)
						continue;
					if (allPlanar || aromaticGeometry(r))
						planarRings.insert(r);
				}
				//~ // find ring systems
				//~ vector<vector<Molib::Fragmenter::Ring>> componentRings = Molib::Fragmenter(residue).identify_fused_rings(planarRings);
				// find ring systems
				//~ set<Molib::Fragmenter::Ring> seenMolib::Fragmenter::Rings;
				//~ vector<Molib::BondSet > fusedMolib::Bonds;
				//~ vector<set<Molib::Atom*> > fusedMolib::Atoms;
				vector<vector<Molib::Ring>> componentRings;
				//~ set<Molib::Atom*> ringAssignedNs;
				for (auto &r : planarRings) {
					bool add_new_cring = true;
					for (auto &cring : componentRings) {
						for (auto &ring : cring) {
							if (ring != r) {
								Molib::Atom::Set inter;
								set_intersection(ring.begin(), ring.end(), r.begin(), r.end(), inserter(inter, inter.begin()));
								if (inter.size() > 1) {
									add_new_cring = false;
									cring.push_back(r);
								}
							}
						}
					}
					if (add_new_cring)
						//~ componentRings.push_back(vector<Molib::Fragmenter::Ring>(r.begin(), r.end()));
						componentRings.push_back(vector<Molib::Ring>{r});
				}
				//~ for (Molib::Fragmenter::Rings::iterator ri = planarRings.begin(); ri != planarRings.end();
											//~ ++ri) {
					//~ const Molib::Fragmenter::Ring &r = *ri;
					//~ if (seenMolib::Fragmenter::Rings.find(r) != seenMolib::Fragmenter::Rings.end())
						//~ continue;
					//~ Molib::BondSet systemMolib::Bonds;
					//~ set<Molib::Atom*> systemMolib::Atoms;
					//~ vector<Molib::Fragmenter::Ring> systemRings;
					//~ vector<Molib::Fragmenter::Ring> queue;
					//~ queue.push_back(r);
					//~ seenMolib::Fragmenter::Rings.insert(r);
					//~ while (queue.size() > 0) {
						//~ Molib::Fragmenter::Ring qr = queue.back();
						//~ queue.pop_back();
						//~ const Molib::Fragmenter::Ring::Molib::Bonds &bonds = qr.bonds();
						//~ const Molib::Fragmenter::Ring::Molib::Atoms &atoms = qr.atoms();
						//~ systemMolib::Bonds.insert(bonds.begin(), bonds.end());
						//~ systemMolib::Atoms.insert(atoms.begin(), atoms.end());
						//~ systemRings.push_back(qr);
						//~ for (Molib::Fragmenter::Ring::Molib::Bonds::const_iterator bi = bonds.begin();
						//~ bi != bonds.end(); ++bi) {
							//~ const Molib::Bond *b = *bi;
							//~ std::vector<Molib::Fragmenter::Ring> bRings = b->minimumMolib::Fragmenter::Rings();
							//~ for (std::vector<Molib::Fragmenter::Ring>::iterator bri =
							//~ bRings.begin(); bri != bRings.end(); ++bri) {
								//~ Molib::Fragmenter::Ring br = *bri;
								//~ if (seenMolib::Fragmenter::Rings.find(br) != seenMolib::Fragmenter::Rings.end())
									//~ continue;
								//~ if (planarRings.find(br) == planarRings.end())
									//~ continue;
								//~ queue.push_back(br);
								//~ seenMolib::Fragmenter::Rings.insert(br);
							//~ }
						//~ }
					//~ }
					//~ fusedMolib::Bonds.push_back(systemMolib::Bonds);
					//~ fusedMolib::Atoms.push_back(systemMolib::Atoms);
					//~ componentRings.push_back(systemRings);
				//~ }
				//~ for (unsigned int i=0; i < fusedMolib::Bonds.size(); ++i) {
				for (auto &systemRings : componentRings) {
					//~ std::set<Molib::Bond *> &bonds = fusedMolib::Bonds[i];
					//~ std::set<Molib::Atom *> &atoms = fusedMolib::Atoms[i];
					Molib::Atom::Set atoms;
					for (auto &ring : systemRings)
						for (auto &a : ring)
							atoms.insert(a);
					Molib::BondSet bonds;
					for (auto &pa : atoms) {
						Molib::Atom &a = *pa;
						//~ for (auto &adj_a : a)
							//~ if (a.atom_number() < adj_a.atom_number())
								//~ bonds.insert(&bond_map.at(Molib::BondKey(&a, &adj_a)));
						//~ for (auto &bond : a) bonds.insert(&bond);
						for (auto &pbond : a.get_bonds()) bonds.insert(pbond);
					}
							//~ else
								//~ bonds.insert(__bond.at({&adjacent, &atom}));
					//~ std::vector<Molib::Fragmenter::Ring> &systemRings = componentRings[i];
			
					if (atoms.size() > 50) {
						// takes too long to do a massive fused-ring system;
						// assume aromatic
						for (auto &pfa : atoms) {
							Molib::Atom &fa = *pfa;
						//~ for (std::set<Molib::Atom *>::iterator fai =
						//~ atoms.begin(); fai != atoms.end(); ++fai) {
							//~ Molib::Atom *fa = *fai;
							if (fa.element() == Molib::Element::C) {
								fa.set_idatm_type("Car");
							} else if (fa.element() == Molib::Element::O) {
								fa.set_idatm_type("Oar");
							}
						}
						continue;
					}
		
					//~ // skip mapped rings
					//~ if (mapped[*atoms.begin()])
						//~ // since rings shouldn't span residues, one atom mapped => all mapped
						//~ continue;
		
					// find bonds directly connected to rings
					// and try to judge their order
					//~ std::map<Molib::Bond *, BondOrder> connected;
					map<Molib::Bond*, BondOrder> connected;
					//~ std::set<std::pair<Molib::Atom *, Molib::Bond*> > ringNeighbors;
					set<pair<Molib::Atom*, Molib::Bond*>> ringNeighbors;
					vector<pair<Molib::Bond*, Molib::Atom*> > possiblyAmbiguous;
					for (auto &pa : atoms) {
						Molib::Atom &a = *pa;
					//~ for (std::set<Molib::Atom *>::iterator ai = atoms.begin();
									//~ ai != atoms.end(); ++ai) {
						//~ Molib::Atom *a = *ai;
						//~ primary = a->primaryNeighbors();
						for (auto &n : a) {
						//~ for (auto &bond : a) {
							//~ auto &n = bond.second_atom();
						//~ for (std::vector<Molib::Atom *>::const_iterator ni =
						//~ primary.begin(); ni != primary.end(); ++ni) {
							//~ Molib::Atom *n = *ni;
							//~ if (atoms.find(n) != atoms.end())
							if (atoms.find(&n) != atoms.end())
								continue;
							//~ Molib::Bond *nb = (*a->bondsMap().find(n)).second;
							//~ Molib::Bond &nb = bond_map.at((a.atom_number() < n.atom_number() ? make_pair(&a, &n)
												//~ : make_pair(&n, &a)));
							//~ Molib::Bond &nb = a[&n];
							Molib::Bond &nb = a.get_bond(n);
							//~ ringNeighbors.insert(std::pair<Molib::Atom *, Molib::Bond *>(n, nb));
							ringNeighbors.insert({&n, &nb});
			
							//~ if (ambiguousVal2Cs.find(n) != ambiguousVal2Cs.end()) {
							if (ambiguousVal2Cs.count(&n)) {
								//~ connected[nb] = AMBIGUOUS;
								connected[&nb] = AMBIGUOUS;
								continue;
							}
							
							//~ Molib::Atom::IdatmInfoMap::const_iterator gi =
								//~ infoMap.find(n->idatmType().str());
							auto gi = help::infoMap.find(n.idatm_type_unmask());
							if (gi == help::infoMap.end()) {
								connected[&nb] = SINGLE;
								continue;
							}
							int geom = (*gi).second.geometry;
							if (geom != 3) {
								connected[&nb] = SINGLE;
								if (geom == 4
								//~ && n->primaryNeighbors().size() == 1) {
								&& n.size() == 1) {
									//~ possiblyAmbiguous.push_back(
										//~ std::pair<Molib::Bond *, Molib::Atom*>(nb, n));
									possiblyAmbiguous.push_back({&nb, &n});
								}
								continue;
							}
							if (n.element() == Molib::Element::N) {
								// aniline can be planar
								connected[&nb] = SINGLE;
								continue;
							}
							// look at neighbors (and grandneighbors)
							bool outsideDouble = false;
							bool ambiguous = (redo[&n] == -1);
							//~ std::vector<Molib::Atom *> nn = n->primaryNeighbors();
							//~ if (nn.size() == 1) {
							if (n.size() == 1) {
								if (a.element() == Molib::Element::N && n.element() == Molib::Element::O)
									connected[&nb] = SINGLE;
								else {
									connected[&nb] = DOUBLE;
									//~ possiblyAmbiguous.push_back(
										//~ std::pair<Molib::Bond *, Molib::Atom*>(nb, n));
									possiblyAmbiguous.push_back({&nb, &n});
								}
								continue;
							}
							//~ for (std::vector<Molib::Atom *>::iterator n2i =
							//~ nn.begin(); n2i != nn.end(); ++n2i) {
								//~ Molib::Atom *n2 = *n2i;
							for (auto &n2 : n) {
							//~ for (auto &bond2 : n) {
								//~ auto &n2 = bond2.second_atom();
								if (&n2 == &a)
									continue;
								auto gi = help::infoMap.find(n2.idatm_type_unmask());
								//~ gi = infoMap.find(
										//~ n2->idatmType().str());
								if (gi == help::infoMap.end())
								//~ if (gi == infoMap.end())
									continue;
								int n2geom = (*gi).second.geometry;
								if (n2geom != 3)
									continue;
								bool allSingle = true;
								//~ std::vector<Molib::Atom *> nnn =
										//~ n2->primaryNeighbors();
								//~ for (std::vector<Molib::Atom *>::iterator n3i
								//~ = nnn.begin(); n3i != nnn.end(); ++n3i){
									//~ Molib::Atom *n3 = *n3i;
								for (auto &n3 : n2) {
								//~ for (auto &bond3 : n2) {
									//~ auto &n3 = bond3.second_atom();
									if (&n3 == &n)
										continue;
									auto gi = help::infoMap.find(n3.idatm_type_unmask());
									//~ gi = infoMap.find(
										//~ n3->idatmType().str());
									//~ if (gi == infoMap.end())
									if (gi == help::infoMap.end())
										continue;
									int n3geom = (*gi).second.geometry;
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
								connected[&nb] = SINGLE;
							else if (ambiguous)
								connected[&nb] = AMBIGUOUS;
							else
								connected[&nb] = DOUBLE;
						}
					}
					//~ std::map<Molib::Bond *, int> curAssign;
					//~ std::vector<std::map<Molib::Bond *, int> > assignments;
					//~ std::vector<std::vector<Molib::Atom *> > assignedUncertains;
					//~ std::map<Molib::Atom *, Molib::Bond *> uncertain2bond;
					map<Molib::Bond*, int> curAssign;
					vector<map<Molib::Bond*, int> > assignments;
					vector<Molib::Atom::Vec > assignedUncertains;
					map<Molib::Atom*, Molib::Bond*> uncertain2bond;
					//~ makeAssignments(bonds, connected, curAssign, &assignments);
					//~ makeAssignments(bonds, connected, curAssign, &assignments,
							//~ bond_map);
					makeAssignments(bonds, connected, curAssign, &assignments);
					if (assignments.size() == 0)
						// try a charged ring
						//~ makeAssignments(bonds, connected, curAssign, &assignments, true);
						//~ makeAssignments(bonds, connected, curAssign, &assignments, 
								//~ bond_map, true);
						makeAssignments(bonds, connected, curAssign, &assignments, true);
					else {
						// if there are no aromatic assignments for a ring and the ring
						// has a nitrogen/oxygen, append charged assignments
						bool addCharged = false;
						//~ for (std::vector<Molib::Fragmenter::Ring>::iterator ri = systemRings.begin();
						//~ ri != systemRings.end(); ++ri) {
							//~ Molib::Fragmenter::Ring ring = *ri;
						for (auto &ring : systemRings) {
							bool hasNO = false;
							//~ for (Molib::Fragmenter::Ring::Molib::Atoms::const_iterator rai = ring.atoms().begin();
							//~ rai != ring.atoms().end(); ++rai) {
								//~ Molib::Atom *a = *rai;
							for (auto &pa : ring) {
								Molib::Atom &a = *pa;
								if (a.element() == Molib::Element::N || a.element() == Molib::Element::O) {
									hasNO = true;
									break;
								}
							}
							if (!hasNO)
								continue;
							bool anyAro = false;
							//~ Molib::BondSet ring_bonds;
							//~ for (auto &pa : ring) {
								//~ Molib::Atom &a = *pa;
								//~ for (auto &adj_a : a)
									//~ if (ring.count(&adj_a))
										//~ if (a.atom_number() < adj_a.atom_number())
											//~ ring_bonds.insert(&bond_map.at({&a, &adj_a}));
							//~ }
							//~ Molib::BondVec rb = get_bonds_in(ring);
							//~ Molib::BondSet ring_bonds(rb.begin(), rb.end());
							Molib::BondSet ring_bonds = get_bonds_in(ring);
							//~ for (std::vector<std::map<Molib::Bond *, int> >::iterator ai =
							//~ assignments.begin(); ai != assignments.end(); ++ai) {
							for (vector<map<Molib::Bond *, int> >::iterator ai =
							assignments.begin(); ai != assignments.end(); ++ai) {
								int bondSum = 0;
								//~ for (Molib::Fragmenter::Ring::Molib::Bonds::const_iterator rbi = ring.bonds().begin();
								//~ rbi != ring.bonds().end(); ++rbi) {
									//~ Molib::Bond *b = *rbi;
								for (auto &b : ring_bonds) {
									//~ bondSum += (*ai)[b];
									bondSum += (*ai)[b];
								}
								//~ int size = ring.bonds().size();
								int size = ring_bonds.size();
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
							//~ std::vector<std::map<Molib::Bond *, int> > prevAssignments = assignments;
							vector<map<Molib::Bond*, int> > prevAssignments = assignments;
							assignments.clear();
							//~ makeAssignments(bonds, connected, curAssign, &assignments, 
									//~ bond_map, true);
							makeAssignments(bonds, connected, curAssign, &assignments, true);
							assignments.insert(assignments.end(), prevAssignments.begin(),
								prevAssignments.end());
						}
					}
					if (assignments.size() == 0) {
						// see if flipping a possibly-ambiguous bond
						// allows an assignment to be made
						//~ std::vector<std::pair<Real, Molib::Bond*> > sortable;
						vector<pair<double, Molib::Bond*> > sortable;
						//~ for (std::vector<std::pair<Molib::Bond *, Molib::Atom*> >::iterator
						//~ si = possiblyAmbiguous.begin();
						//~ si != possiblyAmbiguous.end(); ++si) {
						for (auto &kv : possiblyAmbiguous) {
							//~ Molib::Bond *b = (*si).first;
							//~ Molib::Atom *a = (*si).second;
							Molib::Bond &b = *kv.first;
							Molib::Atom &a = *kv.second;
							Molib::Element e = a.element();
							//~ Real certainty;
							double certainty;
							if (e == Molib::Element::O) {
								//~ certainty = b->sqlength() - p3o2c2;
								//~ certainty = b.first_atom().crd().distance_sq(b.second_atom().crd()) - p3o2c2;
								certainty = b.atom1().crd().distance_sq(b.atom2().crd()) - p3o2c2;
							} else if (e == Molib::Element::S) {
								//~ certainty = b->sqlength() - p3s2c2;
								//~ certainty = b.first_atom().crd().distance_sq(b.second_atom().crd()) - p3s2c2;
								certainty = b.atom1().crd().distance_sq(b.atom2().crd()) - p3s2c2;
							} else {
								certainty = 0.0;
							}
							if (certainty < 0.0)
								certainty = 0.0 - certainty;
							//~ sortable.push_back(std::pair<Real, Molib::Bond *>(
								//~ certainty, b));
							sortable.push_back({certainty, &b});
						}
						//~ std::sort(sortable.begin(), sortable.end());
						sort(sortable.begin(), sortable.end());
						//~ std::vector<Molib::Bond *> flippable;
						Molib::BondVec flippable;
						//~ for (std::vector<std::pair<Real, Molib::Bond*> >::iterator si
						//~ = sortable.begin(); si != sortable.end(); ++si) {
						for (auto &kv : sortable) {
							//~ flippable.push_back((*si).second);
							flippable.push_back(kv.second);
						}
						if (flippable.size() > 0) {
							//~ flipAssign(flippable, atoms, bonds, connected, &assignments);
							//~ flipAssign(flippable, atoms, bonds, connected, &assignments,
										//~ bond_map);
							flipAssign(flippable, atoms, bonds, connected, &assignments);
							if (assignments.size() == 0)
								//~ flipAssign(flippable, atoms, bonds,
									//~ connected, &assignments, true);
								//~ flipAssign(flippable, atoms, bonds,
									//~ connected, &assignments, bond_map, true);
								flipAssign(flippable, atoms, bonds,
									connected, &assignments, true);
						}
					}
		
					if (assignments.size() == 0) {
						// if adjacent carbons were uncertain (i.e. had
						// "redo" values) try changing their type
						//~ std::vector<Molib::Atom *> uncertain;
						Molib::Atom::Vec uncertain;
						//~ for (std::set<std::pair<Molib::Atom *, Molib::Bond *> >::iterator rni =
						//~ ringNeighbors.begin(); rni != ringNeighbors.end(); ++rni) {
						for (auto &kv : ringNeighbors) {
							//~ Molib::Atom *rna = (*rni).first;
							//~ Molib::Bond *rnb = (*rni).second;
							Molib::Atom &rna = *kv.first;
							Molib::Bond &rnb = *kv.second;
							//~ if (redo.find(&rna) != redo.end() && redo[rna] != 0) {
							if (redo.count(&rna) && redo[&rna] != 0) {
								uncertain.push_back(&rna);
								uncertain2bond[&rna] = &rnb;
							}
						}
						if (uncertain.size() > 0) {
							//~ uncertainAssign(uncertain, uncertain2bond, bonds, connected,
										//~ &assignments, &assignedUncertains);
							//~ uncertainAssign(uncertain, uncertain2bond, bonds, connected,
										//~ &assignments, &assignedUncertains, bond_map);
							uncertainAssign(uncertain, uncertain2bond, bonds, connected,
										&assignments, &assignedUncertains);
							if (assignments.size() == 0)
								//~ uncertainAssign(uncertain, uncertain2bond, bonds, connected,
											//~ &assignments, &assignedUncertains, true);
								//~ uncertainAssign(uncertain, uncertain2bond, bonds, connected,
											//~ &assignments, &assignedUncertains, bond_map, true);
								uncertainAssign(uncertain, uncertain2bond, bonds, connected,
											&assignments, &assignedUncertains, true);
						}
					}
					if (assignments.size() == 0) {
						std::cerr << "Cannot find consistent set of"
						" bond orders for ring system containing atom "
						//~ << (*atoms.begin())->oslIdent() << "\n";
						<< (*atoms.begin())->atom_number() << "\n";
						continue;
					}
		
					//~ std::map<Molib::Bond *, int> *bestAssignment = findBestAssignment(
											//~ assignments, systemRings);
					//~ map<Molib::Bond *, int> *bestAssignment = findBestAssignment(
											//~ assignments, systemRings, bond_map);
					map<Molib::Bond *, int> *bestAssignment = findBestAssignment(
											assignments, systemRings);
					if (bestAssignment != NULL && assignedUncertains.size() > 0) {
						//~ unsigned int baIndex = std::find(assignments.begin(),
							//~ assignments.end(), *bestAssignment) - assignments.begin();
						unsigned int baIndex = find(assignments.begin(),
							assignments.end(), *bestAssignment) - assignments.begin();
						invertUncertains(assignedUncertains[baIndex],
										uncertain2bond, &connected);
					}
		
					// see if individual rings are aromatic -- if not
					// then assign types as per best assignment
					//~ for (std::vector<Molib::Fragmenter::Ring>::iterator ri = systemRings.begin();
					//~ ri != systemRings.end(); ++ri) {
						//~ Molib::Fragmenter::Ring ring = *ri;
					for (auto &ring : systemRings) {
						int minFreeElectrons = 0, maxFreeElectrons = 0;
						//~ std::vector<Molib::Bond *> otherFused;
						Molib::BondVec otherFused;
						//~ Molib::BondVec bv = get_bonds_in(ring);
						//~ Molib::BondSet bonds(bv.begin(), bv.end());
						Molib::BondSet bonds = get_bonds_in(ring);
						//~ for (auto &pa : ring) {
							//~ Molib::Atom &a = *pa;
							//~ for (auto &adj_a : a) 
								//~ if (ring.count(&adj_a))
									//~ if (a.atom_number() < adj_a.atom_number())
										//~ bonds.insert(&bond_map.at({&a, &adj_a}));
						//~ }
						//~ const Molib::Fragmenter::Ring::Molib::Bonds &bonds = ring.bonds();
						//~ std::set<Molib::Bond *> ringBonds;
						Molib::BondSet ringBonds;
						ringBonds.insert(bonds.begin(), bonds.end());
						//~ const Molib::Fragmenter::Ring::Molib::Atoms &atoms = ring.atoms();
						bool aro = true;
						//~ for (Molib::Fragmenter::Ring::Molib::Atoms::const_iterator ai = atoms.begin();
						//~ ai != atoms.end(); ++ai) {
							//~ Molib::Atom *a = *ai;
						for (auto &pa : ring) {
							Molib::Atom &a = *pa;
							//~ const Molib::Atom::Molib::Bonds &a_bonds = a->primaryBonds();
							Molib::BondVec a_bonds;
							//~ for (auto &adj_a : a) {
							//~ for (auto &bond : a) {
							for (auto &pbond : a.get_bonds()) {
								//~ Molib::Atom &adj_a = bond.second_atom();
								//~ a_bonds.push_back(&bond_map.at((a.atom_number() < adj_a.atom_number() ?
								//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a))));
								//~ a_bonds.push_back(&bond);
								a_bonds.push_back(pbond);
							}
							//~ Molib::BondSet a_bonds;
							//~ for (auto &adj_a : a)
								//~ a_bonds.insert(bond_map.at({&a, &adj_a}));
							minFreeElectrons++; // a few exceptions below
							if (a_bonds.size() == 2) {
								int element = a.element().number();
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
								Molib::Bond *outBond;
								//~ for (Molib::Atom::Molib::Bonds::const_iterator abi =
								//~ a_bonds.begin(); abi != a_bonds.end();
								//~ ++abi) {
									//~ Molib::Bond *b = *abi;
								for (auto &pb : a_bonds) {
									Molib::Bond &b = *pb;
									if (ringBonds.find(&b)
									!= ringBonds.end())
										continue;
									outBond = &b;
									if (connected.find(&b)
									!= connected.end()) {
										bo = connected[&b];
										break;
									}
									//~ for (std::vector<std::map<Molib::Bond *, int> >::iterator ai = assignments.begin(); ai != assignments.end(); ++ai) {
									for (vector<map<Molib::Bond *, int> >::iterator ai = 
										assignments.begin(); ai != assignments.end(); ++ai) {
										int assignment=(*ai)[&b];
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
								//~ int element = a->element().number();
								int element = a.element().number();
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
										//~ if (outBond->otherMolib::Atom(a)->primaryNeighbors().size() > 1) {
										//~ Molib::Atom &other_atom = (&outBond->first_atom() == &a ? outBond->second_atom()
												//~ : outBond->first_atom());
										Molib::Atom &other_atom = outBond->second_atom(a);
										if (other_atom.size() > 1) {
											bool isFused = false;
											//~ for (std::vector<Molib::Fragmenter::Ring>::iterator sri = systemRings.begin(); sri != systemRings.end(); ++sri) {
											for (auto &sring : systemRings) {
												//~ if (sri == ri) {
												if (sring == ring) {
													continue;
												}
												//~ Molib::Fragmenter::Ring sring = *sri;
												//~ const Molib::Fragmenter::Ring::Molib::Atoms &srMolib::Atoms = sring.atoms();
												//~ if (srMolib::Atoms.find(a) != srMolib::Atoms.end()) {
												if (sring.find(&a) != sring.end()) {
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
						//~ std::set<Molib::Atom *> protonatableNs;
						set<Molib::Atom*> protonatableNs;
						//~ for (Molib::Fragmenter::Ring::Molib::Atoms::const_iterator ai = atoms.begin();
						//~ ai != atoms.end(); ++ai) {
							//~ Molib::Atom *a = *ai;
						for (auto &pa : atoms) {
							Molib::Atom &a = *pa;
							//~ Molib::Atom *a = *ai;
							Molib::BondVec a_bonds;
							//~ for (auto &adj_a : a) {
							//~ for (auto &bond : a) {
							for (auto &pbond : a.get_bonds()) {
								//~ a_bonds.push_back(&bond_map.at((a.atom_number() < adj_a.atom_number() ?
								//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a))));
								//~ a_bonds.push_back(&bond);
								a_bonds.push_back(pbond);
							}
							//~ Molib::BondVec a_bonds;
							//~ for (auto &adj_a : a)
								//~ if (a.atom_number() < adj_a.atom_number())
									//~ a_bonds.push_back(__bond.at({&a, &adj_a}));
							Molib::Element e = a.element();
							if (e == Molib::Element::N) {
								//~ if (a->primaryBonds().size() == 2)
								if (a_bonds.size() == 2)
									protonatableNs.insert(&a);
								else if (bestAssignment != NULL) {
									// distinguish N2+/Npl
									//~ if (aro && isN2plus(bestAssignment, a->primaryBonds())) {
									if (aro && isN2plus(bestAssignment, a_bonds)) {
										// non-ring bond == 0
										//~ a->setComputedIdatmType("N2+");
										a.set_idatm_type("N2+");
										// type bonded isolated O as O3-
										//~ if (a->primaryBonds().size() == 3) {
										if (a_bonds.size() == 3) {
											//~ Molib::Atom::Molib::BondsMap::const_iterator end = a->bondsMap().end();
											//~ for (Molib::Atom::Molib::BondsMap::const_iterator bmi = a->bondsMap().begin();
											//~ bmi != end; ++bmi) {
											for (auto &nb : a) {
											//~ for (auto &bond : a) {
												//~ auto &nb = bond.second_atom();
												//~ Molib::Atom *nb = (*bmi).first;
												if (nb.element() == Molib::Element::O
												//~ && nb->bondsMap().size() == 1)
												&& nb.size() == 1)
													//~ nb->setComputedIdatmType("O3-");
													nb.set_idatm_type("O3-");
											}
										}
									}
								}
								continue;
							}
							if (!aro)
								continue;
							if (e == Molib::Element::C) {
								//~ a->setComputedIdatmType("Car");
								a.set_idatm_type("Car");
							} else if (e == Molib::Element::O) {
								//~ if (isOarPlus(bestAssignment, a->primaryBonds()))
								if (isOarPlus(bestAssignment, a_bonds))
									a.set_idatm_type("Oar+");
								else
									//~ a->setComputedIdatmType("Oar");
									a.set_idatm_type("Oar");
							} else if (e == Molib::Element::S) {
								//~ a->setComputedIdatmType("Sar");
								a.set_idatm_type("Sar");
							} else if (e == Molib::Element::P) {
								//~ if (a->primaryBonds().size() == 2)
								if (a_bonds.size() == 2)
									if (aroEval % 4 != 0)
										aroEval++;
							} else if (e.number() > 20) {
								if (aroEval % 4 != 0)
									aroEval++;
							}
							dbgmsg("pass 4.5 : " << a.idatm_type_unmask());
						}
						if (bestAssignment == NULL && aro) {
							// N2/Npl depends on number of needed electrons
							while (aroEval % 4 != 0 && protonatableNs.size() > 0) {
								aroEval++;
								Molib::Atom *longestN = NULL;
								float bestVal = 0.0;
								for (std::set<Molib::Atom *>::iterator ai =
								protonatableNs.begin(); ai !=
								protonatableNs.end(); ++ai) {
									//~ Molib::Atom *a = *ai;
									Molib::Atom &a = **ai;
									float val = 0.0;
									Molib::BondVec a_bonds;
									//~ for (auto &adj_a : a) {
									//~ for (auto &bond : a) {
									for (auto &pbond : a.get_bonds()) {
										//~ a_bonds.push_back(&bond_map.at((a.atom_number() < adj_a.atom_number() ?
										//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a))));
										//~ a_bonds.push_back(&bond);
										a_bonds.push_back(pbond);
									}
									//~ Molib::BondVec a_bonds;
									//~ for (auto &adj_a : a)
										//~ if (a.atom_number() < adj_a.atom_number())
											//~ a_bonds.push_back(__bond.at({&a, &adj_a}));
									//~ const Molib::Atom::Molib::Bonds &a_bonds = a->primaryBonds();
									//~ for (Molib::Atom::Molib::Bonds::const_iterator bi =
									//~ a_bonds.begin(); bi != a_bonds.end();
									//~ ++bi) {
										//~ Molib::Bond *b = *bi;
									for (auto &b : a_bonds) {
										val += b->length();
										//~ val += b->first_atom().crd().distance(b->second_atom().crd());
										//~ val += b->atom1().crd().distance(b->atom2().crd());
									}
									if (longestN == NULL || val > bestVal) {
										longestN = &a;
										bestVal = val;
									}
									// avoid retyping in pass 7
									redo[&a] = -7;
								}
								//~ longestN->setComputedIdatmType("Npl");
								longestN->set_idatm_type("Npl");
								protonatableNs.erase(longestN);
								ringAssignedNs.insert(longestN);
			
							}
							//~ for (std::set<Molib::Atom *>::iterator ai =
							//~ protonatableNs.begin(); ai != protonatableNs.end();
							//~ ++ai) {
							for (auto &pa : protonatableNs) {
								//~ Molib::Atom *a = *ai;
								Molib::Atom &a = *pa;
								//~ a->setComputedIdatmType("N2");
								a.set_idatm_type("N2");
								ringAssignedNs.insert(&a);
							}
						} else if (bestAssignment != NULL) {
							
							// decide if two-bond nitrogens are N2 or Npl
							for (auto &pa : protonatableNs) {
								Molib::Atom &a = *pa;
							//~ for (std::set<Molib::Atom *>::iterator ai =
							//~ protonatableNs.begin(); ai !=
							//~ protonatableNs.end(); ++ai) {
								//~ Molib::Atom *a = *ai;
								int bondSum = 0;
								Molib::BondVec a_bonds;
								//~ for (auto &adj_a : a) {
								//~ for (auto &bond : a) {
								for (auto &pbond : a.get_bonds()) {
									//~ a_bonds.push_back(&bond_map.at((a.atom_number() < adj_a.atom_number() ?
									//~ make_pair(&a, &adj_a) : make_pair(&adj_a, &a))));
									//~ a_bonds.push_back(&bond);
									a_bonds.push_back(pbond);
								}
								//~ Molib::BondVec a_bonds;
								//~ for (auto &adj_a : a)
									//~ if (a.atom_number() < adj_a.atom_number())
										//~ a_bonds.push_back(__bond.at({&a, &adj_a}));
								//~ const Molib::Atom::Molib::Bonds &a_bonds = a->primaryBonds();
								//~ for (Molib::Atom::Molib::Bonds::const_iterator bi =
								//~ a_bonds.begin(); bi != a_bonds.end();
								//~ ++bi) {
								for (auto &pb : a_bonds) {
									//~ Molib::Bond *b = *bi;
									Molib::Bond &b = *pb;
									//~ bondSum += (*bestAssignment)[b];
									bondSum += (*bestAssignment)[&b];
								}
								if (bondSum == 2)
									//~ a->setComputedIdatmType("Npl");
									a.set_idatm_type("Npl");
								else
									//~ a->setComputedIdatmType("N2");
									a.set_idatm_type("N2");
								ringAssignedNs.insert(&a);
								dbgmsg("pass 4.5 : " << a.idatm_type_unmask());
							}
						}
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
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
	
			//~ if (a->idatmType() != "C2")
			if (a.idatm_type_unmask() != "C2")
				continue;
	
			if (mapped[&a])
				continue;
			
			bool c2possible = false;
			//~ primary = a->primaryNeighbors();
			Molib::Atom::Vec nbValence1;
			//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
			  //~ bi != primary.end(); ++bi) {
				//~ Molib::Atom *bondee = *bi;
			for (auto &bondee : a) {
			//~ for (auto &bond : a) {
				//~ auto &bondee = bond.second_atom();
				//~ Symbol bondeeType = bondee->idatmType();
				const string bondeeType = bondee.idatm_type_unmask();
				//~ Molib::Atom::IdatmInfoMap::const_iterator i = infoMap.find(
								//~ bondeeType.str());
				auto i = help::infoMap.find(bondeeType);
				if (i == help::infoMap.end())
					continue;
				if ((*i).second.geometry == 3 && bondeeType != "Cac"
						&& bondeeType != "N2+"
						// Npl with two bonds or less may be N2
						&& !(bondeeType == "Npl"
						//~ && bondee->primaryNeighbors().size() > 2
						&& bondee.size() > 2
						// because Ng+ isn't assigned until next pass
						//~ && heavys[bondee] > 1)) {
						&& heavys[&bondee] > 1)) {
					c2possible = true;
				  	break;
				//~ } else if (bondee->primaryNeighbors().size() == 1)
				} else if (bondee.size() == 1)
					nbValence1.push_back(&bondee);
			}
	
			if (!c2possible) {
				//~ if (primary.size() == 3 && nbValence1.size() > 0) {
				if (a.size() == 3 && nbValence1.size() > 0) {
					Molib::Atom *best = NULL;
					float bestRatio;
					const char *bestSp2Type;
					//~ for (std::vector<Molib::Atom *>::iterator nbi = nbValence1.begin();
					//~ nbi != nbValence1.end(); ++nbi) {
						//~ Molib::Atom *nb = *nbi;
					for (auto &pnb : nbValence1) {
						Molib::Atom &nb = *pnb;
						const char *sp2Type;
						//~ Real test;
						double test;
						//~ if (nb->element() == Molib::Element::C) {
						if (nb.element() == Molib::Element::C) {
							sp2Type = "C2";
							test = p3c2c;
						} else if (nb.element() == Molib::Element::O) {
							sp2Type = "O2";
							test = p3o2c2;
						} else if (nb.element() == Molib::Element::N) {
							sp2Type = "N2";
							test = p4n2c;
						} else if (nb.element() == Molib::Element::S) {
							sp2Type = "S2";
							test = p3s2c2;
						} else
							continue;
						//~ Real sqlen = sqdistance(nb->coord(), a->coord());
						double sqlen = nb.crd().distance_sq(a.crd());
						//~ Real ratio = sqlen / test;
						double ratio = sqlen / test;
						if (best == NULL || ratio < bestRatio) {
							best = &nb;
							bestRatio = ratio;
							bestSp2Type = sp2Type;
						}
					}
					if (best != NULL)
						//~ best->setComputedIdatmType(bestSp2Type);
						best->set_idatm_type(bestSp2Type);
					else
						//~ a->setComputedIdatmType("C3");
						a.set_idatm_type("C3");
				} else
					//~ a->setComputedIdatmType("C3");
					a.set_idatm_type("C3");
			}
			dbgmsg("pass 5  : " << a.idatm_type_unmask());
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
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
			
			//~ primary = a->primaryNeighbors();
			//~ const Molib::Element &element = a->element();
			const Molib::Element &element = a.element();
			//~ if (element == Molib::Element::N && a->idatmType() != "N3+") {
			if (element == Molib::Element::N && a.idatm_type_unmask() != "N3+") {
				if (mapped[&a])
					continue;
				//~ if (isN3plusOkay(primary))
				if (isN3plusOkay(a))
					//~ a->setComputedIdatmType("N3+");
					a.set_idatm_type("N3+");
				
			//~ } else if (a->idatmType() == "C2") {
			} else if (a.idatm_type_unmask() == "C2") {
				int numNpls = 0;
				//~ for (std::vector<Molib::Atom *>::const_iterator bi =
				  //~ primary.begin(); bi != primary.end(); ++bi) {
					//~ Molib::Atom *bondee = *bi;
				for (auto &bondee : a) {
				//~ for (auto &bond : a) {
					//~ auto &bondee = bond.second_atom();
					//~ if ((bondee->idatmType() == "Npl"
					if ((bondee.idatm_type_unmask() == "Npl"
					  && !mapped[&bondee])
					  || bondee.idatm_type_unmask() == "Ng+")
						// Ng+ possible through template
						// typing
						numNpls++;
				}
	
				bool noplus = false;
				if (numNpls >= 2) {
					//~ for (std::vector<Molib::Atom *>::const_iterator bi =
					  //~ primary.begin(); bi != primary.end(); ++bi) {
						//~ Molib::Atom *bondee = *bi;
					for (auto &bondee : a) {
					//~ for (auto &bond : a) {
						//~ auto &bondee = bond.second_atom();
	
						//~ if (bondee->idatmType() != "Npl")
						if (bondee.idatm_type_unmask() != "Npl")
							continue;
						if (mapped[&bondee])
							continue;
						
						//~ if (bondee->rings(false,
								//~ ringLimit, &considerMapped).size() > 0) {
						Molib::Rings rs = Molib::Fragmenter(residue.get_atoms()).identify_rings();
						Molib::Atom::Set ring_atoms;
						for (auto &ring : rs) 
							for (auto &pa : ring)
								ring_atoms.insert(pa);
						//~ if (bondee.rings(false,
								//~ ringLimit, &considerMapped).size() > 0) {
						if (ring_atoms.count(&bondee)) {
							noplus = true;
							break;
						}
						//~ bondee->setComputedIdatmType("Ng+");
						bondee.set_idatm_type("Ng+");
					}
				}
				if (noplus) {
					for (auto &bondee : a) {
					//~ for (auto &bond : a) {
						//~ auto &bondee = bond.second_atom();
					//~ for (std::vector<Molib::Atom *>::const_iterator bi =
					  //~ primary.begin(); bi != primary.end(); ++bi) {
						//~ Molib::Atom *bondee = *bi;
						if (mapped[&bondee])
							continue;
						if (bondee.idatm_type_unmask() == "Ng+")
							//~ bondee->setComputedIdatmType("Npl");
							bondee.set_idatm_type("Npl");
					}
				}
			//~ } else if (a->idatmType() == "Cac") {
			} else if (a.idatm_type_unmask() == "Cac") {
				//~ for (std::vector<Molib::Atom *>::const_iterator bi =
				  //~ primary.begin(); bi != primary.end(); ++bi) {
					//~ Molib::Atom *bondee = *bi;
				for (auto &bondee : a) {
				//~ for (auto &bond : a) {
					//~ auto &bondee = bond.second_atom();
					if (mapped[&bondee])
						continue;
					if (bondee.element() == Molib::Element::O &&
					  heavys[&bondee] == 1) {
						//~ bondee->setComputedIdatmType("O2-");
						bondee.set_idatm_type("O2-");
					}
				}
			}
			dbgmsg("pass 6  : " << a.idatm_type_unmask());
		}
	
		// "pass 7":  a non-IDATM pass:  split off heavy-atom-valence-2
		//	Npls that have no hydrogens as type N2.
		//	Discrimination criteria is the average bond length of the two 
		//	heavy-atom bonds (shorter implies more double-bond character,
		//	thereby no hydrogen).
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
	
			if (mapped[&a])
				continue;
	
			//~ if (a->idatmType() != "Npl")
			if (a.idatm_type_unmask() != "Npl")
				continue;
			
			if (heavys[&a] != 2)
				continue;
	
			if (ringAssignedNs.find(&a) != ringAssignedNs.end())
				continue;
			
			//~ primary = a->primaryNeighbors();
			//~ if (primary.size() > 2)
			if (a.size() > 2)
				continue;
	
			//~ Real threshold = 1.0;
			//~ Real harmLen = 1.0;
			double threshold = 1.0;
			double harmLen = 1.0;
			Molib::Atom *recipient = NULL;
			//~ Real bratio = 1.0;
			double bratio = 1.0;
			//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
			//~ bi != primary.end(); ++bi) {
				//~ Molib::Atom *bondee = *bi;
			for (auto &bondee : a) {
			//~ for (auto &bond : a) {
				//~ auto &bondee = bond.second_atom();
				//~ Real criteria;
				double criteria;
				if (bondee.element() == Molib::Element::C) {
					criteria = p7cn2nh;
				} else if (bondee.element() == Molib::Element::N) {
					criteria = p7nn2nh;
				} else if (bondee.element() == Molib::Element::O) {
					//~ if (bondee->primaryNeighbors().size() > 1)
					if (bondee.size() > 1)
						continue;
					criteria = p7on2nh;
				} else {
					continue;
				}
				threshold *= criteria;
				//~ Real len = distance(bondee->coord(), a->coord());
				double len = bondee.crd().distance(a.crd());
				harmLen *= len;
				if (len > criteria)
					continue;
				//~ Real ratio = len / criteria;
				double ratio = len / criteria;
				if (ratio > bratio)
					continue;
				//~ if (bondee->element() == Molib::Element::N && bondee->primaryBonds().size() > 2)
				if (bondee.element() == Molib::Element::N && bondee.size() > 2)
					continue;
				if (bondee.idatm_type_unmask() == "Car")
					continue;
	
				//~ std::vector<Molib::Atom *> bondeePrimary = bondee->primaryNeighbors();
				bool doubleOkay = true;
				//~ for (std::vector<Molib::Atom *>::const_iterator bi = bondeePrimary.begin();
				  //~ bi != bondeePrimary.end(); ++bi) {
					//~ Molib::Atom *grandBondee = *bi;
				for (auto &grandBondee : bondee) {
				//~ for (auto &gbond : bondee) {
					//~ auto &grandBondee = gbond.second_atom();
					if (&grandBondee == &a)
						continue;
					//~ Symbol gbType = grandBondee->idatmType();
					const string gbType = grandBondee.idatm_type_unmask();
	
					auto i = help::infoMap.find(gbType);
					//~ Molib::Atom::IdatmInfoMap::const_iterator i = infoMap.find(gbType.str());
					if (i == help::infoMap.end())
					//~ if (i == infoMap.end())
						continue;
					int geom = (*i).second.geometry;
					if (geom > 1 && geom < 4 && heavys[&grandBondee] == 1) {
							doubleOkay = false;
							break;
					}
				}
				if (!doubleOkay)
					continue;
				recipient = &bondee;
				bratio = ratio;
			}
	
			if (harmLen < threshold && recipient != NULL) {
				//~ a->setComputedIdatmType("N2");
				a.set_idatm_type("N2");
				if (recipient->element() == Molib::Element::C) {
					//~ recipient->setComputedIdatmType("C2");
					recipient->set_idatm_type("C2");
				} else if (recipient->element() == Molib::Element::N) {
					//~ recipient->setComputedIdatmType("N2");
					recipient->set_idatm_type("N2");
				} else if (recipient->element() == Molib::Element::O) {
					//~ recipient->setComputedIdatmType("O2");
					recipient->set_idatm_type("O2");
				}
			}
			dbgmsg("pass 7  : " << a.idatm_type_unmask());
		}
	
		// "pass 8":  another non-IDATM: change planar nitrogens bonded only
		//  SP3 atoms to N3 or N3+.  Change Npls/N3s bonded to sp2 atoms that
		//  are not in turn bonded to sp2 atoms (implying Npl doubled bonded)
		//  to N2, otherwise to Npl.  
		//~ std::map<Molib::Atom *, std::vector<Molib::Atom *> > bondedSp2s;
		map<Molib::Atom*, Molib::Atom::Vec> bondedSp2s;
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
	
			if (mapped[&a])
				continue;
	
			if (ringAssignedNs.find(&a) != ringAssignedNs.end())
				continue;
			
			//~ Symbol idatmType = a->idatmType();
			const string idatmType = a.idatm_type_unmask();
			if (idatmType != "Npl" && idatmType != "N2" && idatmType!="N3")
				continue;
			
			bool aroRing = false;
			//~ const Molib::Atom::Molib::Fragmenter::Rings &aRings = a->rings(false, ringLimit, &considerMapped);
			Molib::Rings rs = Molib::Fragmenter(residue.get_atoms()).identify_rings();
			Molib::Rings aRings;
			for (auto &ring : rs) 
				if (ring.count(&a))
					aRings.insert(ring);
			//~ const Molib::Fragmenter::Rings &aRings = a.rings(false, ringLimit, &considerMapped);
			//~ for (Molib::Atom::Molib::Fragmenter::Rings::const_iterator ri = aRings.begin();
			//~ ri != aRings.end(); ++ ri) {
			for (auto &ring : aRings) {
				//~ if ((*ri)->aromatic()) {
				bool aromatic_ring = true; // test that ring is aromatic, i.e., all Carbons are Car
				for (auto &pa : ring) {
					if (pa->element() == Molib::Element::C && pa->idatm_type_unmask() != "Car")
						aromatic_ring = false;
				}
				if (aromatic_ring) {
					aroRing = true;
					break;
				}
			}
			if (aroRing)
				continue;
	
			// any sp2?
			//~ primary = a->primaryNeighbors();
			//~ if (idatmType == "Npl" && primary.size() != 2)
			if (idatmType == "Npl" && a.size() != 2)
				continue;
	
			std::vector<Molib::Atom *> bsp2list;
			//~ for (std::vector<Molib::Atom *>::const_iterator bi = primary.begin();
			//~ bi != primary.end(); ++bi) {
				//~ Molib::Atom *bondee = *bi;
			for (auto &bondee : a) {
			//~ for (auto &bond : a) {
				//~ auto &bondee = bond.second_atom();
				//~ Symbol idatmType = bondee->idatmType();
				const string idatmType = bondee.idatm_type_unmask();
	
				aroRing = false;
				//~ const Molib::Atom::Molib::Fragmenter::Rings &bRings = bondee->rings(false, // ????????????????????????????????
									//~ ringLimit, &considerMapped);
				Molib::Rings rs = Molib::Fragmenter(residue.get_atoms()).identify_rings();
				Molib::Rings bRings;
				for (auto &ring : rs) 
					if (ring.count(&bondee))
						bRings.insert(ring);
				//~ for (Molib::Atom::Molib::Fragmenter::Rings::const_iterator ri = bRings.begin();
				//~ ri != bRings.end(); ++ ri) {
				for (auto &ring : bRings) {
				//~ if ((*ri)->aromatic()) {
					bool aromatic_ring = true; // test that ring is aromatic, i.e., all Carbons are Car
					for (auto &pa : ring) {
						if (pa->element() == Molib::Element::C && pa->idatm_type_unmask() != "Car")
							aromatic_ring = false;
					}
					if (aromatic_ring) {
					//~ if ((*ri)->aromatic()) {
						aroRing = true;
						break;
					}
				}
				if (aroRing) {
					if (heavys[&a] == 1) { // aniline
						//~ a->setComputedIdatmType("Npl");
						a.set_idatm_type("Npl");
						break;
					}
					continue;
				}
	
				auto i = help::infoMap.find(idatmType);
				//~ Molib::Atom::IdatmInfoMap::const_iterator i = infoMap.find(
								//~ idatmType.str());
				if (i == help::infoMap.end() || (*i).second.geometry != 3
						|| bondee.idatm_type_unmask() == "Npl")
					continue;
				bsp2list.push_back(&bondee);
			}
			bondedSp2s[&a] = bsp2list;
		}
	
		// order typing by easiest-figure-out (1 sp2 bonded) to hardest (3)
		// good test cases:  1CY in 3UM8; WRA in 1J3I
		for (unsigned int i=1; i<4; ++i) {
			//~ for (map<Molib::Atom *, vector<Molib::Atom *> >::const_iterator sp2i = bondedSp2s.begin();
			//~ sp2i != bondedSp2s.end(); ++sp2i) {
				//~ const std::vector<Molib::Atom *> &sp2s = (*sp2i).second;
			for (auto &kv : bondedSp2s) {
				const Molib::Atom::Vec &sp2s = kv.second;
				if (sp2s.size() != i)
					continue;
				//~ Molib::Atom *a = (*sp2i).first;
				Molib::Atom &a = *kv.first;
				bool anySP2 = false;
	
				//~ for (std::vector<Molib::Atom *>::const_iterator bi = sp2s.begin();
				//~ bi != sp2s.end(); ++bi) {
					//~ Molib::Atom *bondee = *bi;
				for (auto &pbondee : sp2s) {
					Molib::Atom &bondee = *pbondee;
					anySP2 = true;
					bool remoteSP2 = false;
					//~ std::vector<Molib::Atom *> grandPrimary = bondee->primaryNeighbors();
					//~ for (std::vector<Molib::Atom *>::const_iterator gbi =
					//~ grandPrimary.begin(); gbi != grandPrimary.end(); ++gbi){
					for (auto &grandBondee : bondee) {
					//~ for (auto &gbond : bondee) {
						//~ auto &grandBondee = gbond.second_atom();
						//~ Molib::Atom *grandBondee = *gbi;
						if (&grandBondee == &a)
							continue;
						auto gi = help::infoMap.find(grandBondee.idatm_type_unmask());
						//~ Molib::Atom::IdatmInfoMap::const_iterator gi =
						//~ infoMap.find(grandBondee->idatmType().str());
						if (gi == help::infoMap.end()
						|| (*gi).second.geometry == 3) {
							if (grandBondee.idatm_type_unmask() != "Car"
							&& grandBondee.idatm_type_unmask() != "Npl"
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
						//~ a->setComputedIdatmType("N2");
						a.set_idatm_type("N2");
						break;
					}
					// a remote sp2 atom doesn't necessarily mean Npl
					// (see N1 in MX1 of 2aio), so no else clause
				}
				if (!anySP2) {
					int hvys = heavys[&a];
					if (hvys > 1)
						//~ a->setComputedIdatmType("N3");
						a.set_idatm_type("N3");
					else if (hvys == 1)
						//~ a->setComputedIdatmType(isN3plusOkay(primary) ? "N3+" : "N3");
						a.set_idatm_type(isN3plusOkay(a) ? "N3+" : "N3");
					else
						//~ a->setComputedIdatmType("N3+");
						a.set_idatm_type("N3+");
				}
				dbgmsg("pass 8  : " << a.idatm_type_unmask());
			}
		}
	
		// "pass 9":  another non-IDATM pass and analogous to pass 8:
		//  change O3 bonded only to non-Npl sp2 atom not in turn bonded
		//  to non-Npl sp2 to O2.
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
	
			if (mapped[&a])
				continue;
	
			//~ Symbol idatmType = a->idatmType();
			const string idatmType = a.idatm_type_unmask();
			if (idatmType != "O3")
				continue;
			
			//~ primary = a->primaryNeighbors();
			//~ if (primary.size() != 1)
			if (a.size() != 1)
				continue;
	
			// any sp2?
			//~ Molib::Atom *bondee = *(primary.begin());
			Molib::Atom &bondee = a[0];
			//~ Molib::Atom &bondee = a.first().second_atom();
			//~ Symbol bondeeType = bondee->idatmType();
			const string bondeeType = bondee.idatm_type_unmask();
	
			bool aroRing = false;
			//~ const Molib::Atom::Molib::Fragmenter::Rings &bRings = bondee->rings(false, ringLimit, &considerMapped);			
			Molib::Rings rs = Molib::Fragmenter(residue.get_atoms()).identify_rings();
			Molib::Rings bRings;
			for (auto &ring : rs) 
				if (ring.count(&bondee))
					bRings.insert(ring);
			//~ for (Molib::Atom::Molib::Fragmenter::Rings::const_iterator ri = bRings.begin();
			//~ ri != bRings.end(); ++ ri) {
			for (auto &ring : bRings) {
				bool aromatic_ring = true; // test that ring is aromatic, i.e., all Carbons are Car
				for (auto &pa : ring) {
					if (pa->element() == Molib::Element::C && pa->idatm_type_unmask() != "Car")
						aromatic_ring = false;
				}
				if (aromatic_ring) {
				//~ if (ri.aromatic()) {
					aroRing = true;
					break;
					//~ aroRing = true;
					//~ break;
				}
			}
				//~ }
			//~ }
			if (aroRing) {
				// can't be O2
				continue;
			}
	
			auto i = help::infoMap.find(bondeeType);
			//~ Molib::Atom::IdatmInfoMap::const_iterator i = infoMap.find(
								//~ bondeeType.str());
			if (i == help::infoMap.end() || (*i).second.geometry != 3)
				continue;
			bool remoteSP2 = false;
			//~ std::vector<Molib::Atom *> grandPrimary = bondee->primaryNeighbors();
			//~ for (std::vector<Molib::Atom *>::const_iterator gbi =
			//~ grandPrimary.begin(); gbi != grandPrimary.end(); ++gbi) {
				//~ Molib::Atom *grandBondee = *gbi;
			for (auto &grandBondee : bondee) {
			//~ for (auto &gbond : bondee) {
				//~ auto &grandBondee = gbond.second_atom();
				if (&grandBondee == &a)
					continue;
				auto gi = help::infoMap.find(grandBondee.idatm_type_unmask());
				//~ Molib::Atom::IdatmInfoMap::const_iterator gi =
					//~ infoMap.find(grandBondee->idatmType().str());
				if (gi == help::infoMap.end() || (*gi).second.geometry == 3) {
					if (grandBondee.idatm_type_unmask() != "Car"
					&& grandBondee.idatm_type_unmask() != "Npl") {
						remoteSP2 = true;
						break;
					}
				}
			}
			if (!remoteSP2)
				//~ a->setComputedIdatmType("O2");
				a.set_idatm_type("O2");
			dbgmsg("pass 9  : " << a.idatm_type_unmask());
		}
	
		// "pass 10":  another non-IDATM pass. Ensure nitrate ions are N2/O2-
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain)
		for (auto &a : residue) {
		//~ for (Molib::Atoms::iterator ai = pMolib::Atoms.begin(); ai != pMolib::Atoms.end(); ++ai) {
			//~ Molib::Atom *a = *ai;
	
			if (mapped[&a])
				continue;
	
			if (a.element() != Molib::Element::N)
				continue;
			
			//~ primary = a->primaryNeighbors();
			//~ if (primary.size() != 2)
			if (a.size() != 2)
				continue;
	
			bool bondersOkay = true;
			//~ for (Molib::Atoms::const_iterator pi = primary.begin(); pi != primary.end(); ++pi) {
				//~ Molib::Atom *bondee = *pi;
			for (auto &bondee : a) {	
			//~ for (auto &bond : a) {	
				//~ auto &bondee = bond.second_atom();
				if (bondee.element() != Molib::Element::O
				//~ || bondee->primaryNeighbors().size() != 1) {
				|| bondee.size() != 1) {
					bondersOkay = false;
					break;
				}
			}
	
			if (bondersOkay) {
				//~ a->setComputedIdatmType("N2");
				a.set_idatm_type("N2");
				for (auto &bondee : a) {
				//~ for (auto &bond : a) {
					//~ auto &bondee = bond.second_atom();
				//~ for (Molib::Atoms::const_iterator pi = primary.begin(); pi != primary.end(); ++pi)
					//~ (*pi)->setComputedIdatmType("O2-");
					bondee.set_idatm_type("O2-");
				}
			}
			dbgmsg("pass 10  : " << a.idatm_type_unmask());
		}
		dbgmsg(molecule);
	}
};
