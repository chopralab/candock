#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "molecule.hpp"
#include "geom3d/geom3d.hpp"
#include "fragmenter/fragmenter.hpp"
#include "fragmenter/unique.hpp"
#include "openmm/forcefield.hpp"
#include "atomtype.hpp"
using namespace std;

namespace Molib {
	int freeOxygens(const Atom &a, map<const Atom*, int> &heavys) {
		int freeOxygens = 0;
		for (auto &bondee : a) {
			if (bondee.element() == Element::O && heavys[&bondee] == 1)
				freeOxygens++;
		}
		return freeOxygens;
	}
	
	bool aromatic(const Atom::Set &ring) {
		double sum = 0.0;
		int bonds = 0;
		dbgmsg("ring = " << endl << ring);
		auto ring_bonds = get_bonds_in(ring);
		dbgmsg("aromaticity test for ring bonds : " << endl 
			<< ring_bonds << endl << "--------------------");
		for (auto &pbond : ring_bonds) {
			Bond &bond = *pbond;
			const Atom &atom1 = bond.atom1();
			const Atom &atom2 = bond.atom2();
			const Element &e1 = atom1.element();
			const Element &e2 = atom2.element();
			double d = bond.length();
			dbgmsg(e1 << " " << e2);
			if (e1 == Element::C && e2 == Element::C) {
				bonds++;
				sum += (d - 1.397) * (d - 1.397);
			} else if ((e1 == Element::C || e2 == Element::C) &&
			  (e1 == Element::N || e2 == Element::N)) {
				bonds++;
				sum += (d - 1.338) * (d - 1.338);
			} else if (e1 == Element::N && e2 == Element::N) {
				bonds++;
				sum += (d - 1.308) * (d - 1.308);
			} else if ((e1 == Element::C || e2 == Element::C) &&
			  (e1 == Element::O || e2 == Element::O)) {
				bonds++;
				sum += (d - 1.300) * (d - 1.300);
			}
		}
		if (bonds == 0)
			return false;
		double homas = 1.0 - (98.89 / bonds) * sum;
		dbgmsg(homas);
		if (homas < 0.271) {
			// not aromatic
			return false;
		}
		return true;
	}
	void AtomType::compute_idatm_type(Molecule &molecule) {
		// angle values used to discriminate between hybridization states
		const double angle23val1 = 115.0;
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
		const double p3o2c2 = 1.30 * 1.30;
		const double p3o2as = 1.685 * 1.685;
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
	
		const double p8cc2n2h = 1.367 * 1.367;
		const double p8nn2n2h = 1.326 * 1.326;
		const double p8cn2n2h = 1.367 * 1.326;
	
		// algorithm based on E.C. Meng / R.A. Lewis paper 
		// "Determination of Molecular Topology and Atomic Hybridization
		// States from Heavy Atom Coordinates", J. Comp. Chem., v12#7, 891-898
		// and on example code from idatm.f implementation by E.C. Meng
	
		// differences: No boron types.  Double-bonded Npls are split off
		//   as N2.  Sox split into Sxd (sulfoxide), and Son (sulfone).
		//   Carbons in aromatic rings are type Car.  Aromatic oxygens are Oar.
	
		// initialize idatm type in Atoms
		for (auto &pa : molecule.get_atoms())
			pa->set_idatm_type("???");
	
		map<const Atom*, int> heavys; // number of heavy atoms bonded
	
		// "pass 1":  type hydrogens / deuteriums and compute number of
		// heavy atoms connected to each atom
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			const Element &element = a.element();
			dbgmsg(a.element());
			if (element.number() == 1) {
				// sort out if it's a hydrogen or deuterium
				bool isHyd = true;
				if (a.atom_name()[0] == 'D' || a.atom_name()[0] == 'd')
					isHyd = false;
				bool bondedToCarbon = false;
				for (auto &bondee : a) {
					if (bondee.element() == Element::C) {
						bondedToCarbon = true;
						break;
					}
				}
				a.set_idatm_type(bondedToCarbon ? (isHyd ? "HC" : "DC") : (isHyd ? "H" : "D"));
				dbgmsg("hydrogens" << boolalpha << (bondedToCarbon ? (isHyd ? "HC" : "DC") : (isHyd ? "H" : "D")));
			}
	
			int heavyCount = 0;
			for (auto &bondee : a) {
				if (bondee.element().number() > 1) {
					heavyCount++;
				}
			}
			heavys[&a] = heavyCount;
			dbgmsg("pass 1  : " << a.idatm_type_unmask() << " heavys[" << a.atom_name() << "] = " << heavys[&a]);
		}
	
		// "pass 1.5": use templates for "infallible" typing of standard residue types
		map<const Atom*, bool> mapped;
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model) {
			for (auto &residue : chain) {
				auto it = help::standard_residues.find(residue.resn());
				if (it != help::standard_residues.end())
					for (auto &a : residue) {
						if (!mapped[&a]) {
							string idatm_type = "???";
							// is it the N-terminal residue ?
							if (residue.rest() == Residue::protein 
								&& &residue == &chain.first()) {
								if (a.atom_name() == "N")
									idatm_type = "N3+";
								else if (a.atom_name() == "H1" || a.atom_name() == "H2" || 
										a.atom_name() == "H3" || a.atom_name() == "HN1" || 
										a.atom_name() == "HN2" || a.atom_name() == "HN3")
									idatm_type = "H";
							}
							// is it C-terminal ?
							if (residue.rest() == Residue::protein
								&& &residue == &chain.last()) {
								if (a.atom_name() == "O" || a.atom_name() == "OXT")
									idatm_type = "O2-";
								else if (a.atom_name() == "C") 
									idatm_type = "Cac";
							}
							// H on peptide N
							if (a.atom_name() == "HN")
								idatm_type = "H";
							// still not known? find it in the table!
							if (idatm_type == "???") {
								auto it2 = it->second.find(a.atom_name());
								if (it2 == it->second.end())
									throw Error("die : cannot find atom name " 
										+ a.atom_name() + " of template residue "
										+ residue.resn());
								idatm_type = it2->second;
							}
							a.set_idatm_type(idatm_type);
							mapped[&a] = true;
						}
						dbgmsg("pass 1.5  : " << a.idatm_type_unmask() << " " << a.atom_name() << " " << a.atom_number());
				}
			}
		}
		// "pass 2": elements that are typed only by element type
		// and valences > 1
		map<Atom*, int> redo;
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			if (mapped[&a])
				continue;
			const Element &element = a.element();
	
			// undifferentiated types
			if (element >= Element::He && element <= Element::Be
			  || element >= Element::Ne && element <= Element::Si
			  || element >= Element::Cl) {
			  	a.set_idatm_type(element.name());
			  	dbgmsg(a.idatm_type_unmask());
				continue;
			}
	
			// valence 4
			//	C must be sp3 (C3)
			//	N must be part of an N-oxide (Nox) or a quaternary
			//		amine (N3+) 
			//	P must be part of a phosphate (Pac), a P-oxide (Pox)
			//		or a quaternary phosphine (P3+)
			//	S must be part of a sulfate, sulfonate or sulfamate
			//		(Sac), or sulfone (Son)
			if (a.size() == 4) {
				if (element == Element::C) {
					a.set_idatm_type("C3");
				} else if (element == Element::N) {
					a.set_idatm_type(freeOxygens(a, heavys) > 1 ? "Nox" : "N3+");
				} else if (element == Element::P) {
					int freeOxys = freeOxygens(a, heavys);
					if (freeOxys >= 2)
						a.set_idatm_type("Pac");
					else if (freeOxys == 1)
						a.set_idatm_type("Pox");
					else
						a.set_idatm_type("P3+");
				} else if (element == Element::S) {
					int freeOxys = freeOxygens(a, heavys);
					if (freeOxys >= 3) {
						a.set_idatm_type("Sac");
					} else if (freeOxys >= 1) {
						a.set_idatm_type("Son");
					} else {
						a.set_idatm_type("S");
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
				avgAngle /= 3.0;
	
				if (element == Element::C) {
					if (avgAngle < angle23val1)
						a.set_idatm_type("C3");
					else
						a.set_idatm_type(freeOxygens(a, heavys) >= 2 ?  "Cac" : "C2");
				} else if (element == Element::N) {
					if (avgAngle < angle23val1)
						a.set_idatm_type("N3");
					else
						a.set_idatm_type(freeOxygens(a, heavys) >= 2 ?  "Ntr" : "Npl");
				} else if (element == Element::S) {
					bool hasOxy = false;
					for (int i = 0; i < 3; ++i) {
						if (a[i].element() == Element::O) {
						  	hasOxy = true;
							break;
						}
					}
					a.set_idatm_type(hasOxy ? "Sxd" : "S3+");
				}
			}
	
			// valence 2
			// calculate the bond angle and assign a tentative atom
			// type accordingly (a single angle is often not a good
			// indicator of type).  Mark these atoms for further
			// analysis by putting a non-zero value for them in the
			// 'redo' array.
			//	C may be sp3 (C3), sp2 (C2), or sp (C1)
			//	N may be sp3 (N3), sp2 or planar (Npl), or sp (N1)
			//	O and S are sp3 (O3 and S3, respectively)
			else if (a.size() == 2) {
				double ang = Geom3D::degrees(Geom3D::angle(a[0].crd(), a.crd(), a[1].crd()));
	
				if (element == Element::C) {
					if (ang < angle23val1) {
						a.set_idatm_type("C3");
						redo[&a] = 1;
					} else if (ang < angle12val) {
						a.set_idatm_type("C2");
						if (ang < angle23val2) {
							redo[&a] = 3;
						}
					} else {
						a.set_idatm_type("C1");
					}
				} else if (element == Element::N) {
					if (ang < angle23val1) {
						a.set_idatm_type("N3");
						redo[&a] = 2;
					} else {
						a.set_idatm_type(ang < angle12val ? "Npl" : "N1");
					}
				} else if (element == Element::O) {
					a.set_idatm_type("O3");
				} else if (element == Element::S) {
					a.set_idatm_type("S3");
				}
			}
	
			if (a.idatm_type_unmask() == "???") {
				a.set_idatm_type(element.name());
			}
			dbgmsg("pass 2  : " << a.idatm_type_unmask());
		}
	
		// "pass 3": determine types of valence 1 atoms.  These were typed
		// by element only in previous pass, but can be typed more accurately
		// now that the atoms they are bonded to have been typed.  Bond
		// lengths are used in this pass.  
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			if (a.size() != 1)
				continue;
			
			Atom &bondee = a[0];
			double sqlen = bondee.crd().distance_sq(a.crd());
			string bondeeType = bondee.idatm_type_unmask();
			
			if (a.idatm_type_unmask() == "C") {
				if (mapped[&a])
					continue;
				if (sqlen <= p3c1c1 && bondeeType == "C1") {
					a.set_idatm_type("C1");
				} else if (sqlen <= p3c2c &&
				  bondee.element() == Element::C) {
					a.set_idatm_type("C2");
				} else if (sqlen <= p3c2n &&
				  bondee.element() == Element::N) {
					a.set_idatm_type("C2");
				} else {
					a.set_idatm_type("C3");
				}
			} else if (a.idatm_type_unmask() == "N") {
				if (mapped[&a])
					continue;
				if (sqlen <= p3n1c1 && bondeeType == "C1") {
					a.set_idatm_type("N1");
				} else if (sqlen > p3n3c &&
				  (bondeeType == "C2" || bondeeType == "C3")) {
					a.set_idatm_type("N3");
				} else if ((sqlen > p3n3n3 && bondeeType == "N3") || (sqlen > p3n3n2 && bondeeType == "Npl")) {
					a.set_idatm_type("N3");
				} else {
					a.set_idatm_type("Npl");
				}
			} else if (a.idatm_type_unmask() == "O") {
				if (bondeeType == "Cac" || bondeeType == "Pac" || bondeeType == "Sac" || bondeeType == "Ntr") {
					if (!mapped[&a])
						//~ a.set_idatm_type("O-");
						a.set_idatm_type("O2-");
				} else if (bondeeType == "Nox" || bondeeType == "Pox" || bondeeType == "Son" || bondeeType == "Sxd") {
					if (!mapped[&a])
						a.set_idatm_type("O2");
				} else if (sqlen <= p3o2c2 && bondee.element() == Element::C) {
					if (!mapped[&a])
						a.set_idatm_type("O2");
					if (!mapped[&bondee])
						bondee.set_idatm_type("C2");
					redo[&bondee] = 0;
				} else if (sqlen <= p3o2as && bondee.element() == Element::As) {
					if (!mapped[&a])
						a.set_idatm_type("O2");
				} else {
					if (!mapped[&a])
						a.set_idatm_type("O3");
				}
			} else if (a.idatm_type_unmask() == "S") {
				if (bondee.element() == Element::P) {
					if (!mapped[&a])
						a.set_idatm_type("S2");
				} else if (sqlen <= p3s2c2 && bondee.element() == Element::C) {
					if (!mapped[&a])
						a.set_idatm_type("S2");
					if (!mapped[&bondee])
						bondee.set_idatm_type("C2");
					redo[&bondee] = 0;
				} else if (sqlen <= p3s2as &&
				  bondee.element() == Element::As) {
					if (!mapped[&a])
						a.set_idatm_type("S2");
				} else {
					if (!mapped[&a])
						a.set_idatm_type("S3");
				}
			}
			dbgmsg("pass 3  : " << a.idatm_type_unmask());
		}
	
		// "pass 4": re-examine all atoms with non-zero 'redo' values and
		//   retype them if necessary
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			if (mapped[&a])
				redo[&a] = 0;
	
			if (redo[&a] == 0)
				continue;
			
			bool c3able = false;
			for (auto &bondee : a) {
				const Element &bondeeElement = bondee.element();
				double sqlen = bondee.crd().distance_sq(a.crd());
	
				if (redo[&a] == 1) {
					if ((sqlen > p4c3c && bondeeElement == Element::C) || (sqlen > p4c3n && bondeeElement == Element::N) || (sqlen > p4c3o && bondeeElement == Element::O)) {
						a.set_idatm_type("C3");
						break;
					}
					if ((sqlen <= p4c2c && bondeeElement == Element::C) || (sqlen <= p4c2n && bondeeElement == Element::N)) {
						a.set_idatm_type("C2");
					}
				} else if (redo[&a] == 2) {
					if ((sqlen <= p4n2c && bondeeElement == Element::C) || (sqlen <= p4n2n && bondeeElement == Element::N)) {
						a.set_idatm_type("Npl");
						break;
					}
				} else {
					if ((sqlen <= p4c2c && bondeeElement == Element::C) || (sqlen <= p4c2n && bondeeElement == Element::N)) {
						a.set_idatm_type("C2");
						c3able = false;
						break;
					}
					if ((sqlen > p4c3c && bondeeElement == Element::C) || (sqlen > p4c3n && bondeeElement == Element::N) || (sqlen > p4c3o && bondeeElement == Element::O)) {
						c3able = true;
					}
					if (sqlen > p4ccnd && bondeeElement == Element::C) {
						c3able = true;
					}
				}
			}
			if (c3able)
				a.set_idatm_type("C3");
			dbgmsg("pass 4  : " << a.idatm_type_unmask());
		}
	
		// "pass 5": change isolated sp2 carbons to sp3 since it is 
		//   impossible for an atom to be sp2 hybrizided if all its 
		//   neighbors are sp3 hybridized.  In addition, a carbon atom cannot
		//   be double bonded to a carboxylate carbon, phosphate phosphorus,
		//   sulfate sulfur, sulfone sulfur, sulfoxide sulfur, or sp1 carbon.
		//   Addition not in original idatm: Nox also
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			if (a.idatm_type_unmask() != "C2")
				continue;
	
			if (mapped[&a])
				continue;
			
			bool c2possible = false;
			for (auto &bondee : a) {
				const string bondeeType = bondee.idatm_type_unmask();
	
				if (bondeeType != "C3" && bondeeType != "DC" &&
				  bondeeType != "HC" && bondeeType != "N3" &&
				  bondeeType != "N3+" && bondeeType != "O3" &&
				  bondeeType != "Cac" && bondeeType != "Pac" &&
				  bondeeType != "Sac" && bondeeType != "Son" &&
				  bondeeType != "C1" && bondeeType != "S3" &&
				  bondeeType != "Nox" && bondeeType != "Sxd") {
					c2possible = true;
				  	break;
				}
			}
	
			if (!c2possible)
				a.set_idatm_type("C3");
			dbgmsg("pass 5  : " << a.idatm_type_unmask());
		}
	
		// "pass 6": 
		//   1) make decisions about the charge states of nitrogens.  If an
		//      sp3 nitrogen is bonded to sp3 carbons and/or hydrogens and/or
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
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			if (a.idatm_type_unmask() == "N3") {
				if (mapped[&a])
					continue;
				bool positive = true;
				for (auto &bondee : a) {
					const string bondeeType = bondee.idatm_type_unmask();
	
					if (bondeeType != "C3" && bondeeType != "H" && bondeeType != "D") {
					  	positive = false;
						break;
					}
				}
				if (positive)
					a.set_idatm_type("N3+");
				
			} else if (a.idatm_type_unmask() == "C2") {
				int numNpls = 0;
				for (auto &bondee : a) {
					if ((bondee.idatm_type_unmask() == "Npl" && !mapped[&bondee]) || bondee.idatm_type_unmask() == "Ng+")
						// Ng+ possible through template typing
						numNpls++;
				}
	
				bool noplus = false;
				if (numNpls == 3) {
					for (auto &bondee : a) {
	
						if (bondee.idatm_type_unmask() != "Npl")
							continue;
						if (mapped[&bondee])
							continue;
						
						bondee.set_idatm_type("Ng+");
						for (auto &bondee2 : bondee) {

							const string bondee2type = bondee2.idatm_type_unmask();
							if ((bondee2type == "C2" || bondee2type == "Npl") && &bondee2 != &a) {
								bondee.set_idatm_type("Npl");
								noplus = true;
								break;
							}
						}
	
					}
				}
				if (noplus) {
					for (auto &bondee : a) {
						if (mapped[&bondee])
							continue;
						if (bondee.idatm_type_unmask() == "Ng+")
							bondee.set_idatm_type("Npl");
					}
				}
			} else if (a.idatm_type_unmask() == "Cac") {
				for (auto &bondee : a) {
					if (mapped[&bondee])
						continue;
					if (bondee.element() == Element::O && heavys[&bondee] == 1) {
						//~ bondee.set_idatm_type("O-");
						bondee.set_idatm_type("O2-");
					}
				}
			}
			dbgmsg("pass 6  : " << a.idatm_type_unmask());
		}
	
		// "pass 7":  this pass is not in the IDATM paper but is a suggested
		//    improvement mentioned on page 897 of the paper:  find aromatic
		//    ring types.  The method is to:
		//
		//	1) Find all intraresidue rings
		//	2) Check that all the atoms of the ring are planar types
		//	3) Check bond lengths around the ring; see if they are
		//		consistent with aromatic bond lengths
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if (! help::standard_residues.count(residue.resn())) {
				dbgmsg("here I am");
				Fragmenter frag(residue.get_atoms());
				Rings rings = frag.identify_rings();
				for (auto &r : rings) {
					bool planarTypes = true;
					for (auto &pa : r) {
						const Atom &a = *pa;
						string idatmType = a.idatm_type_unmask();
						if (idatmType != "C2" && idatmType != "Npl" &&
						  idatmType != "S2" && idatmType != "O3" &&
						  idatmType != "S3" && idatmType != "N3" &&
						  idatmType != "Oar" && idatmType != "P" &&
						  idatmType != "Car" &&
						  !(idatmType == "C3" && a.size() == 2)) {
							// O3/S3/N3 will be changed to O2/S2/Npl if
							// ring is found to be aromatic
							dbgmsg(idatmType);
							planarTypes = false;
							break;
						}
			
						if ((idatmType == "C2" || idatmType == "C3"
						  || idatmType == "O3" || idatmType == "S3"
						  || idatmType == "N3") && mapped[&a]) {
							dbgmsg(idatmType);
						  	planarTypes = false;
							break;
						}
					}
					dbgmsg("before planarTypes = " << boolalpha << planarTypes);
					if (!planarTypes)
						continue;

					dbgmsg("before aromatic");
					if (!aromatic(r)) {
						dbgmsg("not aromatic");
						continue;
					}

					// Correct types to be aromatic
					for (auto &pa : r) {
						Atom &a = *pa;
						a.add_property("ag" + help::to_string(r.size()))
							.add_property("ag");
						if (a.idatm_type_unmask() == "C2")
							a.set_idatm_type("Car");
						if (a.idatm_type_unmask() == "C3")
							a.set_idatm_type("Car");
						if (a.idatm_type_unmask() == "O3")
							a.set_idatm_type("Oar");
						if (a.idatm_type_unmask() == "S3")
							a.set_idatm_type("S2");
						if (a.idatm_type_unmask() == "N3")
							a.set_idatm_type("Npl");
						dbgmsg("pass 7  : " << a.idatm_type_unmask() << " " << pa->idatm_type());
					}		
				}
			}
		}
		// "pass 8":  another non-IDATM pass:  split off heavy-atom-valence-2
		//	Npls that have no hydrogens as type N2.
		//	Discrimination criteria is the average bond length of the two 
		//	heavy-atom bonds (shorter implies more double-bond character,
		//	thereby no hydrogen).
		for (auto &pa : molecule.get_atoms()) {
			Atom &a = *pa;
			if (mapped[&a])
				continue;
	
			if (a.idatm_type_unmask() != "Npl")
				continue;
			
			if (heavys[&a] != 2)
				continue;
			
			if (a.size() > 2)
				continue;
	
			// are both bonded heavy atoms sp3?  If so -> Npl
			bool bothSP3 = true;
			bool bothC = true;
			bool bothN = true;
			bool other = false;
			for (auto &bondee : a) {
				const string idatmType = bondee.idatm_type_unmask();
	
				if (idatmType != "C3"
				&& idatmType != "N3+"
				&& idatmType != "N3"
				&& idatmType != "O3"
				&& idatmType != "S3+"
				&& idatmType != "S3"
				&& idatmType != "P3+") {
					bothSP3 = false;
					// don't break out, need to check both
					// bonded heavies' element types
				}
				
				if (bondee.element() == Element::C) {
					bothN = false;
				} else if (bondee.element() == Element::N) {
					bothC = false;
				} else {
					other = true;
					bothN = false;
					bothC = false;
				}
			}
			if (bothSP3)
				continue;
	
			if (other) {
				// beats me; guess single bonds
				continue;
			}
			double avgLen = 0.0;
			for (auto &bondee : a) {
				double sqlen = bondee.crd().distance_sq(a.crd());
				avgLen += sqlen;
			}
			avgLen /= 2.0;
			if (bothC) {
				if (avgLen <= p8cc2n2h)
					a.set_idatm_type("N2");
			} else if (bothN) {
				if (avgLen <= p8nn2n2h)
					a.set_idatm_type("N2");
			} else if (avgLen <= p8cn2n2h) {
				// one N, one C
				a.set_idatm_type("N2");
			}
			dbgmsg("pass 8  : " << a.idatm_type_unmask());
		}
		// "pass 9": change O2- to O3- for sulfates, phosphates, N-oxide, S3- for thiophosphate ... 
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if ( ! help::standard_residues.count(residue.resn())) {
				Fragmenter(residue.get_atoms()).substitute_atoms(help::special);
				dbgmsg("pass 9 (renaming special atom types) : " << endl << residue);
			}
		}
		//~ // missing types: C1-,N2+,N1+,Oar+,O1+,O1,Sar
		// still missing types: C1-,O1+,O1
		// no longer missing : Sar,Oar+,N1+,N2+
		dbgmsg("MOLECULE AFTER IDATM TYPING" << endl << molecule);
	}	
	void AtomType::refine_idatm_type(Molecule &molecule) {
		// "pass 10": refined idatm typing relying on bond orders... 
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if ( ! help::standard_residues.count(residue.resn())) {
				Fragmenter(residue.get_atoms()).substitute_atoms(help::refine);
				dbgmsg("pass 10 (refining idatm atom types) : " << endl << residue);
			}
		}
		dbgmsg("MOLECULE AFTER REFINED IDATM TYPING" << endl << molecule);
	}
	
	void AtomType::compute_gaff_type(Molecule &molecule) {
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if (!help::standard_residues.count(residue.resn())) {
				Fragmenter f(residue.get_atoms());
				f.substitute_atoms(help::gaff);
				dbgmsg("pass 1 (gaff applying basic rules) : " << endl << residue);
				f.flip_conjugated_gaff_types(residue.get_atoms());
				dbgmsg("pass 2 (gaff distincting conjugated types "
					<< "cc(cd), ce(cf), nc(nd), ne(nf), pe(pf), and "
					<< "bridge type cp(cq)) : " << endl << residue);
			}
		}
		dbgmsg("MOLECULE AFTER GAFF TYPING" << endl << molecule);
	}

	BondVec get_double_bonds(const BondSet &bonds) {
		BondVec dbl;
		for (auto &pbond : bonds)
			if (pbond->is_double()) dbl.push_back(pbond);
		return dbl;
	}

	bool two_continuous_single_bonds(const BondSet &bonds) {
		for (auto i = bonds.begin(); i != bonds.end(); ++i) {
			if ((*i)->is_single()) {
				auto j = i;
				for (++j; j != bonds.end(); ++j) {
					if ((*j)->is_single()) {
						if ((*i)->is_adjacent(**j)) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}
	
	void AtomType::compute_ring_type(Molecule &molecule) {
		dbgmsg("Computing ring types for molecule " << molecule.name());
		for (auto &assembly : molecule)
		for (auto &model : assembly)
		for (auto &chain : model)
		for (auto &residue : chain) {
			if (! help::standard_residues.count(residue.resn())) {
				for (auto &atom : residue)
					atom.insert_property("NG", 1); // atom belongs to chain
				Fragmenter frag(residue.get_atoms());
				Rings rings = frag.identify_rings();
				//~ Rings rings = frag.identify_fused_rings();
				for (auto &ring : rings) {
					for (auto &pa : ring) {
						pa->erase_property("NG"); // atom belongs to ring
					}
				}
				int r_id = 0;
				for (auto &ring : rings) {
					auto ring_bonds = get_bonds_in(ring);
					auto out_bonds = get_bonds_in(ring, false);
					for (auto &pbond : ring_bonds) {
						Bond &bond = *pbond;
						bond.set_ring(true); // ring bond
					}
					// Pure aliphatic atom in a ring, which is made of sp3 carbon
					string aromatic_type = "AR5"; 
					// determine the type of aromatic ring
					if (aromatic(ring)) {
						int num_c = 0, num_n = 0;
						for (auto &pa : ring) {
							if (pa->idatm_type_unmask() == "Car") num_c++;
							else if (pa->idatm_type_unmask() == "N2" 
								|| pa->idatm_type_unmask() == "N2+") num_n++;
						}
						// calculate number of out-of-ring double bonds
						auto dbl = get_double_bonds(out_bonds);
						int num_double_out = 0;
						for (auto &pbond : dbl) {
							if (pbond->atom1().has_property("NG")
							 || pbond->atom2().has_property("NG"))
								++num_double_out;
						}			 
						dbgmsg("num_c = " << num_c << " num_n = " << num_n
							<< " num_double_out = " << num_double_out);
						if (ring.size() == 6 && num_c + num_n == 6 && num_double_out == 0) {
							// Pure aromatic atom (such as benzene and pyridine)
							aromatic_type = "AR1";
						} else if (num_double_out > 0) {
							// Atom in a planar ring, which has one or several double bonds 
							// formed between non-ring atoms and the ring atoms
							aromatic_type = "AR3";
						} else if (get_double_bonds(ring_bonds).size() >= 2 
							&& two_continuous_single_bonds(ring_bonds)) {
							// Atom in a planar ring, usually the ring has two 
							// continous single bonds and at least two double bonds
							aromatic_type = "AR2";
						} else {
							// Atom other than AR1, AR2, AR3 and AR5.
							aromatic_type = "AR4";
						}
					}
					const string ring_type = "RG" + help::to_string(ring.size());
					for (auto &pa : ring) {
						pa->add_property(ring_type)
							.add_property(aromatic_type);
					}
				}
			}
		}
		dbgmsg("MOLECULE AFTER RING TYPING" << endl << molecule);
	}
};
