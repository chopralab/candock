#ifndef BOND_H
#define BOND_H
#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "it.hpp"
#include "graph/graph.hpp"

using namespace std;

namespace Molib {	
	class Atom;
	//~ typedef pair<Atom*, Atom*> AtomPair;
	
	class Bond : public template_vector_container<Bond*, Bond> {
		Atom *__atom1, *__atom2;
		string __rotatable;
		string __bond_gaff_type;
		int __bo;
		bool __ring;
		bool __owns_atoms;
		set<int> __angles;
		int __drive_id;
	public:
		Bond() : __atom1(nullptr), __atom2(nullptr), __owns_atoms(false),
			__ring(false), __bo(0), __drive_id(0), __rotatable(""),
			__bond_gaff_type("") {}
		//~ Bond(Atom *atom1, Atom *atom2, bool owns_atoms=false) : __atom1(atom1), __atom2(atom2), 
			//~ __owns_atoms(owns_atoms), __ring(false) {}
		Bond(Atom *atom1, Atom *atom2, bool owns_atoms=false) : __atom1(atom1), 
			__atom2(atom2), __owns_atoms(owns_atoms), __ring(false), __bo(0), 
			__drive_id(0), __rotatable(""), __bond_gaff_type("") {}
		//~ Bond(Atom *atom1, Atom *atom2, const string &rotatable) : __atom1(atom1), 
			//~ __atom2(atom2), __owns_atoms(false), __rotatable(rotatable), __ring(false) {}
		//~ Bond(Atom *atom1, Atom *atom2, const string &order) 
			//~ : __atom1(atom1), __atom2(atom2), __bond_gaff_type(order), __owns_atoms(true),
			//~ __ring(false) {}
		//~ ~Bond() { if (__owns_atoms) { if (__atom1) { delete __atom1; __atom1 = nullptr; }
			//~ if (__atom2) { delete __atom2; __atom2 = nullptr; }	} }
		~Bond();
		void set_members(const string &str);
		//~ bool is_adjacent(const Bond &other) { 
			//~ return &first_atom() == &other.first_atom() 
				//~ || &second_atom() == &other.second_atom() 
				//~ || &first_atom() == &other.second_atom() 
				//~ || &second_atom() == &other.first_atom();
		//~ }
		//~ bool is_adjacent(const Bond &other) { 
			//~ return &atom1() == &other.atom1() 
				//~ || &atom2() == &other.atom2() 
				//~ || &atom1() == &other.atom2() 
				//~ || &atom2() == &other.atom1();
		//~ }
		bool is_adjacent(const Bond &other);
		//~ bool is_terminal() const { return __atom1->get_num_heavy() == 1 || __atom2->get_num_heavy() == 1; }
		void set_angle(int angle) { __angles.insert(angle); }
		void set_drive_id(int drive_id) { __drive_id = drive_id; }
		void set_rotatable(const string &rotatable) { __rotatable = rotatable; }
		const string& get_rotatable() const { return __rotatable; }
		void set_ring(bool ring) { __ring = ring; }
		bool is_ring() const { return __ring; }
		bool is_rotatable() const { return !__rotatable.empty(); }
		//~ void set_bond_atoms(Atom &atom1, Atom &atom2) { __atom1 = &atom1; __atom2 = &atom2; }
		//~ void set_bond_type(const string &bond_type) { __bond_type = bond_type; }
		//~ const string& get_bond_type() const { return __bond_type; }
		void set_bond_gaff_type(const string &bond_gaff_type) { __bond_gaff_type = bond_gaff_type; }
		const string& get_bond_gaff_type() const { return __bond_gaff_type; }
		void set_bo(const int &bo) { __bo = bo; }
		int get_bo() const { return __bo; }
		bool is_single() const { return __bo == 1; }
		bool is_double() const { return __bo == 2; }
		bool is_triple() const { return __bo == 3; }
		bool is_aromatic() const { return __bond_gaff_type == "AB" 
			|| __bond_gaff_type == "SB" || __bond_gaff_type == "DB"; }
		//~ Bond& get_reverse() const { return __atom2[&*__atom1]; }
		//~ Bond& get_reverse() const;
		//~ Atom& first_atom() const { return *__atom1; }
		//~ Atom& second_atom() const { return *__atom2; }
		Atom& first_atom(Atom &origin) const { return (&origin == __atom1 ? *__atom1 : *__atom2); }
		Atom& second_atom(Atom &origin) const { return (&origin == __atom1 ? *__atom2 : *__atom1); }
		Atom& atom1() const { return *__atom1; }
		Atom& atom2() const { return *__atom2; }
		//~ double length() const { return first_atom().crd().distance(second_atom().crd()); }
		double length() const;
		//~ bool operator==(const Bond& rhs) { return &first_atom() == &rhs.first_atom() 
			//~ && &second_atom() == &rhs.second_atom() || &first_atom() == &rhs.second_atom() 
			//~ && &second_atom() == &rhs.first_atom(); }
		//~ bool operator!=(const Bond& rhs) { return !(*this == rhs); }
		friend ostream& operator<< (ostream& stream, const Bond& b);
		// the following are required for BondGraph :-)
		bool compatible(const Bond &other) const;
		//~ bool compatible(const Atom &atom) const;
		//~ string get_label() const { return get_bond_gaff_type(); }
		//~ string get_label() const { stringstream ss; ss << *this; return ss.str(); }
		string get_label() const;
		//~ string print() const { return get_label(); }
		const int weight() const { return 0; } // dummy for graph ostream operator
	};
	
	//~ typedef Glib::Graph<Bond> MolGraph;
	
	typedef Glib::Graph<Bond> BondGraph;

	typedef vector<Bond*> BondVec;
	typedef set<Bond*> BondSet;
	map<string, int> decode_smiles_prop(const vector<string> &s);
	vector<unique_ptr<Bond>> create_bonds(const help::smiles &edges);
	//~ MolGraph create_graph(const help::smiles &edges);
	BondGraph create_graph(const help::smiles &edges);
	//~ void connect_bonds(const BondVec &bonds);
	void connect_bonds(const BondSet &bonds);
	//~ void connect_bonds(BondVec &bonds);

	//~ template<typename T>
	//~ MolGraph create_graph(const T &atoms) {
		//~ return MolGraph(get_bonds_in(atoms), true);
	//~ }
	//~ template<typename T>
	//~ BondGraph create_graph(const T &atoms) {
		//~ return BondGraph(get_bonds_in(atoms), true);
	//~ }
	//~ BondGraph create_graph(const BondVec &bonds);
	BondGraph create_graph(const BondSet &bonds);
	ostream& operator<< (ostream& stream, const BondVec& bonds);
	ostream& operator<< (ostream& stream, const BondSet& bonds);
} // Molib
#endif
