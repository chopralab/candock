#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"

namespace Molib {
	
	bool Bond::is_adjacent(const Bond &other) { 
		return atom1().atom_number() == other.atom1().atom_number() 
			|| atom2().atom_number() == other.atom2().atom_number() 
			|| atom1().atom_number() == other.atom2().atom_number() 
			|| atom2().atom_number() == other.atom1().atom_number();
	}

	string Bond::get_label() const { 
		stringstream ss; 
		ss << atom1().get_label() << "#" << atom1().atom_number() << "_" 
			<< atom2().get_label() << "#" << atom2().atom_number()
			<< "_" << get_bond_gaff_type();
		return ss.str();
	}

	Bond::~Bond() { 
		if (__owns_atoms) {
			dbgmsg("deleting bond");
			if (__atom1) { 
				delete __atom1; 
				__atom1 = nullptr; 
			}
			if (__atom2) { 
				delete __atom2; 
				__atom2 = nullptr; 
			}
		}
	}

	double Bond::length() const { 
		return atom1().crd().distance(atom2().crd()); 
	}

	bool Bond::compatible(const Bond &other) const { 
		/* this is "smiles" bond and "other" is real molecule bond */
		return (atom1().compatible(other.atom1()) && atom2().compatible(other.atom2())
			|| atom1().compatible(other.atom2()) && atom2().compatible(other.atom1()))
			&& (get_bond_gaff_type().empty() || get_bond_gaff_type() == other.get_bond_gaff_type())
			&& (get_bo() == 0 || get_bo() == other.get_bo());
	}

	void Bond::set_members(const string &str) {
		// set angles
		boost::match_results<string::const_iterator> matches;
		auto i = str.find("angles");
		if (i != string::npos) {
			string::const_iterator begin = str.begin() + i, end = str.end();
			while (boost::regex_search(begin, end, matches, boost::regex("([-]*\\d+)[,}]"))) {
				string s(matches[1].first, matches[1].second);
				int angle = stoi(s);
				set_angle(angle);
				begin = matches[1].second;
				if (*begin == '}') break;
			}
		}
		// set rotatable type
		boost::smatch m;
		if (boost::regex_search(str, m, boost::regex("rota=([^,$]+)"))) {
			if (m[1].matched) {
				set_rotatable(m[1].str());
			}
		}
		// set drive_id
		if (boost::regex_search(str, m, boost::regex("drive_id=(\\w+)"))) {
			if (m[1].matched) {
				set_drive_id(stoi(m[1].str()));
			}
		}
		// set gaff bond type
		if (boost::regex_search(str, m, boost::regex("bond_gaff_type=(\\w+)"))) {
			if (m[1].matched) {
				set_bond_gaff_type(m[1].str());
			}
		}
		// set bond order
		if (boost::regex_search(str, m, boost::regex("bo=(\\w+)"))) {
			if (m[1].matched) {
				set_bo(stoi(m[1].str()));
			}
		}
	}

	void connect_bonds(const BondSet &bonds) {
		for (auto &pbond1 : bonds) {
			dbgmsg("bond is " << *pbond1);
			BondSet connected;
			for (auto &pbond2 : pbond1->atom1().get_bonds()) {
				dbgmsg("    bond of atom1 = " << *pbond2);
				if (pbond1 != pbond2) connected.insert(pbond2);
			}
			for (auto &pbond2 : pbond1->atom2().get_bonds()) {
				dbgmsg("    bond of atom2 = " << *pbond2);
				if (pbond1 != pbond2) connected.insert(pbond2);
			}
			for (auto &pbond : connected) pbond1->add(pbond);
			dbgmsg("conn = " << endl << connected << endl << "end conn");
		}
		dbgmsg("connected bonds : " << endl << bonds
			<< endl << "----------------");
	}
	
	void erase_bonds(const BondSet &bonds) {
		for (auto &pbond : bonds) {
			pbond->clear();
		}
	}
	
	BondGraph create_graph(const BondSet &bonds) {
		return BondGraph(bonds, true);
	}
	
	map<string, int> decode_smiles_prop(const vector<string> &s) {
		map<string, int> smiles_prop;
		if (s.size() > 2) {
			auto bonds = help::ssplit(s.at(2), ",", true);
			if (bonds.size() > 0) {
				smiles_prop["B"] = stoi(bonds.at(0));
			}
			if (bonds.size() == 2) {
				smiles_prop["H"] = stoi(bonds.at(1).substr(0, 1)); // e.g. 3H
			}
		}
		if (s.size() > 3) { 
			auto vec = help::ssplit(s.at(3), ",", true); 
			for (auto &prop : vec) {
				int cnt = isdigit(prop.at(0)) ? stoi(prop.substr(0, 1)) : -1;
				string nm = isdigit(prop.at(0)) ? prop.substr(1) : prop;
				smiles_prop[nm] = cnt;
			}
		}
		return smiles_prop;
	}

	vector<unique_ptr<Bond>> create_bonds(const help::smiles &edges) {
		map<int, int> added;
		Atom::Vec atoms;
		vector<unique_ptr<Bond>> bonds;
		dbgmsg("creating bond graph from smiles");
		for (auto &e : edges) {
			auto s1 = help::ssplit(e.atom_property1, "#", true);
			auto s2 = help::ssplit(e.atom_property2, "#", true);
			string smiles_label1 = s1.at(0);
			string smiles_label2 = s2.at(0);
			int idx1 = stoi(s1.at(1));
			int idx2 = stoi(s2.at(1));
			if (!added.count(idx1)) {
				added[idx1] = atoms.size();
				map<string, int> smiles_prop1 = decode_smiles_prop(s1);
				atoms.push_back(new Atom(added[idx1] + 1, smiles_label1, smiles_prop1));
				dbgmsg("idx1 = " << idx1 << " added[idx1] = " << added[idx1]);
			}
			if (!added.count(idx2)) {
				added[idx2] = atoms.size();
				map<string, int> smiles_prop2 = decode_smiles_prop(s2);
				atoms.push_back(new Atom(added[idx2] + 1, smiles_label2, smiles_prop2));
				dbgmsg("idx2 = " << idx2 << " added[idx2] = " << added[idx2]);
			}
			Atom &atom1 = *atoms[added[idx1]];
			Atom &atom2 = *atoms[added[idx2]];
			// here bond "owns" atoms and must delete them in destructor
			bonds.push_back(unique_ptr<Bond>(new Bond(new Atom(atom1), new Atom(atom2), true)));
			bonds.back()->set_members(e.bond_property);
		}
		for (int i = 0; i < bonds.size(); ++i) {
			Bond &bond1 = *bonds[i];
			for (int j = i + 1; j < bonds.size(); ++j) {
				Bond &bond2 = *bonds[j];
				if (bond1.is_adjacent(bond2)) {
					bond1.add(&bond2);
					bond2.add(&bond1);
				}
			}
		}
		// must delete pointers after work...
		for (auto &patom : atoms) delete patom;
		return bonds;
	}
	BondGraph create_graph(const help::smiles &edges) {
		return BondGraph(std::move(create_bonds(edges)), true, false);
	}

	void Bond::erase_stale_refs(const Bond &deleted_bond, const BondVec &bonds) {
		for (auto &pbond : bonds) {
			auto &bond = *pbond;
#ifndef NDEBUG
			dbgmsg("bonds before deletion");
			for (auto &adj : bond) dbgmsg(adj);
#endif
			for (int j = 0; j < bond.size(); ++j) {
				if (&deleted_bond == &bond[j]) {
					dbgmsg("deleting reference to bond"
						<< endl << " bond " << bond 
						<< " is no longer connected to bond "
						<< bond[j]);
					bond.erase(j--);
				}
			}
#ifndef NDEBUG
			dbgmsg("bonds after deletion");
			for (auto &adj : bond) dbgmsg(adj);
#endif
		}
	}

	ostream& operator<< (ostream& stream, const Bond& b) {
		return stream << "(address = " << &b << " atom1 = " 
			<< (b.__atom1 ? help::to_string(b.atom1().atom_number()) : "nullptr")
			<< ", atom2 = " << (b.__atom2 ? help::to_string(b.atom2().atom_number()) : "nullptr") 
			<< ", bond_gaff_type = " << b.get_bond_gaff_type() << ", rotatable = " 
			<< b.get_rotatable() << ", aromatic = " << boolalpha << b.is_aromatic() << ")";
	}
	ostream& operator<< (ostream& stream, const BondVec& bonds) {
		for (auto &pbond : bonds) {
			stream << " bond : " << *pbond << endl;
			for (auto &bond2 : *pbond) 
				stream << "       is connected to : " << bond2 << endl;
		}
		return stream;
	}
	ostream& operator<< (ostream& stream, const BondSet& bonds) {
		BondVec bv(bonds.begin(), bonds.end());
		return stream << bv;
	}
}; // Molib
