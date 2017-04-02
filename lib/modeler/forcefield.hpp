#ifndef FORCEFIELD_H
#define FORCEFIELD_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/coordinate.hpp"
#include "helper/debug.hpp"
using namespace std;

namespace Molib {
	class Atom;
	class Molecule;
	class Molecules;
	class Score;
}
namespace OMMIface {
	class ParameterError : public Error {
	public: 
		ParameterError(const string &msg) : Error(msg) {}
	};
	struct ForceField {

		double coulomb14scale;
		double lj14scale;
		double step; // for knowledge-based potential
		
		//~ static const double Temperature         = 300;     // Kelvins
		//~ static const double FrictionInPerPs     = 91.;     // collisions per picosecond
		//~ static const double SolventDielectric   = 80.;     // typical for water
		//~ static const double SoluteDielectric    = 2.;      // typical for protein

		struct AtomType { string cl, element; double mass, 
			charge, sigma, epsilon;	};
		struct BondType { double length, k; bool can_constrain; };
		struct AngleType { double angle, k; };
		struct TorsionType { int periodicity; double phase, k; };
		typedef vector<TorsionType> TorsionTypeVec;
		struct ResidueTopology {
			map<string, int> atom; // name, type
			map<const string, map<const string, bool>> bond; // from, to
			set<string> external_bond; // from
		};
		struct KBType { vector<double> potential, derivative; };
		typedef map<const int, AtomType> Atoms;
		typedef map<const string, map<const string, BondType>> Bonds;
		typedef map<const string, map<const string, map<const string, AngleType>>> Angles;
		typedef map<const string, map<const string, map<const string, map<const string, TorsionTypeVec>>>> Torsions;
		typedef map<const string, ResidueTopology> Residues;
		typedef map<const int, map<const int, KBType>> KBForces;
		Atoms atom_type;
		Bonds bond_type;
		Angles angle_type;
		Torsions torsion_type, improper_type;
		Residues residue_topology;
		KBForces kb_force_type;
                bool has_kb_force_type ( const string &aclass1, const string &aclass2) const;
                bool has_atom_type(      const int aclass1) const;
                bool has_bond_type(      const string &aclass1, const string &aclass2) const;
                bool has_angle_type(     const string &aclass1, const string &aclass2, const string &aclass3) const;
                bool has_dihedral_type(  const string &aclass1, const string &aclass2, const string &aclass3, const string &aclass4) const;
                bool has_improper_type(  const string &aclass1, const string &aclass2, const string &aclass3, const string &aclass4) const;
                const KBType& get_kb_force_type(const Molib::Atom &atom1, const Molib::Atom &atom2) const;
                const AtomType& get_atom_type(const int type) const;
                const BondType& get_bond_type(const int type1, const int type2) const;
                const AngleType& get_angle_type(const int type1, const int type2, const int type3) const;
                const TorsionTypeVec& get_dihedral_type(const int type1, const int type2, const int type3, const int type4) const;
                const TorsionTypeVec& get_improper_type(const int type1, const int type2, const int type3, const int type4) const;
		ForceField& parse_forcefield_file(const string&);
		void output_forcefield_file(const string&);
		ForceField& parse_gaff_dat_file(const string&);
		ForceField& insert_topology(const Molib::Molecule&);
		ForceField& erase_topology(const Molib::Molecule&);
		ForceField& add_kb_forcefield(const Molib::Score&, const double step);
		
		bool residue_exists(const string &name) const { return residue_topology.count(name); }
	};
	ostream& operator<< (ostream& stream, const ForceField::ResidueTopology& r);
	ostream& operator<< (ostream& stream, const ForceField& ff);
}
#endif
