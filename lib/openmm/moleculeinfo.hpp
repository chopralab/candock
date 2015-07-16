#ifndef MOLECULEINFO_H
#define MOLECULEINFO_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/coordinate.hpp"
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "pdbreader/molecule.hpp"
#include "OpenMM.h"
using namespace std;

namespace Molib {
	class Atom;
	class Molecule;
};

namespace OMMIface {
	class ForceField;

	struct MoleculeInfo {
		typedef vector<pair<Molib::Atom*, Molib::Atom*>> Bonds;
		typedef vector<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*>> Angles;
		typedef vector<tuple<Molib::Atom*, Molib::Atom*, Molib::Atom*, Molib::Atom*>> Dihedrals;
		Molib::Atom::Vec atom;
		Bonds bond;
		Angles angle;
		//~ Dihedrals dihedral;
		Dihedrals dihedral, improper;
		Bonds kbforce;
		map<const Molib::Atom*, const int> atom_to_type;
		void __bonds_angles_dihedrals(const Molib::Molecule&);
		~MoleculeInfo() {
			dbgmsg("calling destructor of MoleculeInfo");
		}
		MoleculeInfo& get_molecule_info(const Molib::Molecule&, 
				const ForceField&);
		//~ MoleculeInfo& get_kb_force_info(const Molib::Molecule &, const Molib::Molecule&,
										//~ const ForceField&, const int&);
		MoleculeInfo& get_kb_force_info(const Molib::Molecule &, const Molib::Molecule&,
										const int&);
		MoleculeInfo& get_interaction_force_info(const Molib::Molecule &receptor, 
			const Molib::Molecule &ligand, const int &dist_cutoff);
	};
	ostream& operator<< (ostream& stream, const MoleculeInfo::Bonds &bonds);
	ostream& operator<< (ostream& stream, const MoleculeInfo::Angles &angles);
	ostream& operator<< (ostream& stream, const MoleculeInfo::Dihedrals &dihedrals);

};
#endif
