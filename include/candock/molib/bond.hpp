#ifndef BOND_H
#define BOND_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "candock/molib/it.hpp"
#include "candock/graph/graph.hpp"

using namespace std;

namespace candock {

namespace molib {
        class Atom;
        class Bond;

	class Bond : public template_vector_container<Bond*, Bond> {
		Atom *__atom1, *__atom2;
		int __idx1, __idx2;
	std::string __rotatable;
	std::string __bond_gaff_type;
        std::string __stereo;
		int __bo;
		bool __ring;
		bool __owns_atoms;
		set<int> __angles;
		int __drive_id;
	public:
		Bond() : __atom1(nullptr), __atom2(nullptr), __idx1(0), __idx2(0), __rotatable(""), __bond_gaff_type(""),
			 __bo(0), __ring(false), __owns_atoms(false),  __angles(set<int>()), __drive_id(0) {}
		Bond(Atom *atom1, Atom *atom2, bool owns_atoms=false) : __atom1(atom1), 
			__atom2(atom2), __idx1(0), __idx2(0), __rotatable(""), __bond_gaff_type(""), __bo(0), __ring(false),
                        __owns_atoms(owns_atoms),__angles(set<int>()), __drive_id(0){}
		Bond(Atom *atom1, Atom *atom2, int idx1, int idx2, bool owns_atoms=false) : __atom1(atom1), 
			__atom2(atom2), __idx1(idx1), __idx2(idx2), __rotatable(""), __bond_gaff_type(""), __bo(0), __ring(false),
			 __owns_atoms(owns_atoms), __angles(set<int>()), __drive_id(0) {}
		~Bond();
		bool is_set() const { return __atom1 != nullptr; }
		void set_members(const std::string &str);
		bool is_adjacent(const Bond &other);
		void set_angle(int angle) { __angles.insert(angle); }
		void set_drive_id(int drive_id) { __drive_id = drive_id; }
		void set_rotatable(const std::string &rotatable) { __rotatable = rotatable; }
		const std::string& get_rotatable() const { return __rotatable; }
		void set_ring(bool ring) { __ring = ring; }
		bool is_ring() const { return __ring; }
		bool is_rotatable() const { return !__rotatable.empty() && __rotatable != "amide"; }
		void set_bond_gaff_type(const std::string &bond_gaff_type) { __bond_gaff_type = bond_gaff_type; }
		const std::string& get_bond_gaff_type() const { return __bond_gaff_type; }
		void set_stereo(std::string s) { __stereo = s; }
                const std::string& stereo() const { return __stereo; }
		void set_bo(const int &bo) { __bo = bo; }
		int get_bo() const { return __bo; }
		bool is_single() const { return __bo == 1; }
		bool is_double() const { return __bo == 2; }
		bool is_triple() const { return __bo == 3; }
		bool is_aromatic() const { return __bond_gaff_type == "AB" 
			|| __bond_gaff_type == "SB" || __bond_gaff_type == "DB"; }
		Atom& first_atom(const Atom &origin) const { return (&origin == __atom1 ? *__atom1 : *__atom2); }
		Atom& second_atom(const Atom &origin) const { return (&origin == __atom1 ? *__atom2 : *__atom1); }
		Atom& atom1() const { return *__atom1; }
		Atom& atom2() const { return *__atom2; }
		int idx1() const { return __idx1; }
		int idx2() const { return __idx2; }
		double length() const;
		friend ostream& operator<< (ostream& stream, const Bond& b);
		// the following are required for BondGraph :-)
		bool compatible(const Bond &other) const;
                std::string get_label() const;
                std::string get_pdb_remark() const;
		int weight() const { return 0; } // dummy for graph ostream operator
	};

        struct BondPtrComp
        {
                bool operator()(const Bond* lhs, const Bond* rhs) const;
        };

        typedef std::vector<Bond*> BondVec;
        typedef std::set<Bond*, BondPtrComp> BondSet;
        typedef graph::Graph<Bond> BondGraph;

        void erase_stale_bond_refs(const Bond &deleted_bond, const BondVec &bonds);
	std::map<std::string, int> decode_smiles_prop(const std::vector<std::string> &s);
	std::vector<std::unique_ptr<Bond>> create_bonds(const help::smiles &edges);
	BondGraph create_graph(const help::smiles &edges);
	void connect_bonds(const BondSet &bonds);
	void erase_bonds(const BondSet &bonds);
	BondGraph create_graph(const BondSet &bonds);

	std::ostream& operator<< (std::ostream& stream, const BondVec& bonds);
	std::ostream& operator<< (std::ostream& stream, const BondSet& bonds);
} // molib

}

#endif
