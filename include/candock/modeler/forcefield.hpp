#ifndef FORCEFIELD_H
#define FORCEFIELD_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "candock/geometry/coordinate.hpp"

namespace candock {

namespace molib {
	class Atom;
	class Molecule;
	class Molecules;
}

namespace score {
        class KBFF;
}

namespace OMMIface {
	class ParameterError : public Error {
	public: 
		ParameterError(const std::string &msg) : Error(msg) {}
	};
	struct ForceField {

		double coulomb14scale;
		double lj14scale;
		double step; // for knowledge-based potential
		double kb_cutoff; // also for kb potential

        std::map<const std::string, const int> gaff_name_to_type;

		struct AtomType { std::string cl, element; double mass, 
			charge, sigma, epsilon;	};
		struct BondType { double length, k; bool can_constrain; };
		struct AngleType { double angle, k; };
		struct TorsionType { int periodicity; double phase, k; };
		typedef std::vector<TorsionType> TorsionTypeVec;
		struct ResidueTopology {
			std::map<std::string, int> atom; // name, type
			std::map<const std::string, std::map<const std::string, bool>> bond; // from, to
			std::set<std::string> external_bond; // from
		};
		struct KBType { std::vector<double> potential, derivative; };
		typedef std::map<const int, AtomType> Atoms;
		typedef std::map<const std::string, std::map<const std::string, BondType>> Bonds;
		typedef std::map<const std::string, std::map<const std::string, std::map<const std::string, AngleType>>> Angles;
		typedef std::map<const std::string, std::map<const std::string, std::map<const std::string, std::map<const std::string, TorsionTypeVec>>>> Torsions;
		typedef std::map<const std::string, ResidueTopology> Residues;
		typedef std::map<const int, std::map<const int, KBType>> KBForces;
		Atoms atom_type;
		Bonds bond_type;
		Angles angle_type;
		Torsions torsion_type, improper_type;
		Residues residue_topology;
		KBForces kb_force_type;
                bool has_kb_force_type ( const std::string &aclass1, const std::string &aclass2) const;
                bool has_atom_type(      const int aclass1) const;
                bool has_bond_type(      const std::string &aclass1, const std::string &aclass2) const;
                bool has_angle_type(     const std::string &aclass1, const std::string &aclass2, const std::string &aclass3) const;
                bool has_dihedral_type(  const std::string &aclass1, const std::string &aclass2, const std::string &aclass3, const std::string &aclass4) const;
                bool has_improper_type(  const std::string &aclass1, const std::string &aclass2, const std::string &aclass3, const std::string &aclass4) const;
                const KBType& get_kb_force_type(const molib::Atom &atom1, const molib::Atom &atom2) const;
                const KBType& get_kb_force_type(int aclass1, int aclass2) const;
                const AtomType& get_atom_type(const int type) const;
                const BondType& get_bond_type(const int type1, const int type2) const;
                const AngleType& get_angle_type(const int type1, const int type2, const int type3) const;
                const TorsionTypeVec& get_dihedral_type(const int type1, const int type2, const int type3, const int type4) const;
                const TorsionTypeVec& get_improper_type(const int type1, const int type2, const int type3, const int type4) const;
		ForceField& parse_forcefield_file(const std::string&);
		void output_forcefield_file(const std::string&);
		ForceField& parse_gaff_dat_file(const std::string&);
		ForceField& insert_topology(const molib::Molecule&);
		ForceField& erase_topology(const molib::Molecule&);
		ForceField& add_kb_forcefield(const score::KBFF&, double);
		
		bool residue_exists(const std::string &name) const { return residue_topology.count(name); }
	};
	std::ostream& operator<< (std::ostream& stream, const ForceField::ResidueTopology& r);
	std::ostream& operator<< (std::ostream& stream, const ForceField& ff);
}

}

#endif
