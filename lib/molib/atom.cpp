#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <boost/regex.hpp>
#include "candock/molib/bond.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/residue.hpp"
#include "candock/molib/chain.hpp"
using namespace std;

namespace candock {
namespace molib {
	ostream& operator<< (ostream& stream, const Atom& a) {
		if (!a.__br) {
			stream << "ATOM " << a.get_label() << " " << a.atom_number() << endl;
			return stream;
		}
		stream << setw(6) << left << (a.br().rest() == molib::Residue::hetero 
			|| a.br().rest() == molib::Residue::ion 
			|| a.br().rest() == molib::Residue::water ? "HETATM" : "ATOM");
		stream << setw(5) << right << a.atom_number();
		stream << setw(1) << " ";
		stream << setw(4) << left << (a.atom_name().size() < 4 ? " " + a.atom_name() : a.atom_name().substr(0,4));
		stream << setw(1) << " ";
		stream << setw(3) << right << a.br().resn();
		stream << setw(1) << " ";
		stream << setw(1) << a.br().br().chain_id();
		stream << setw(4) << right << a.br().resi();
		stream << setw(1) << a.br().ins_code();
		stream << setw(27) << a.crd().pdb();
		stream << setw(6) << setprecision(2) << fixed << right << 1.0;
		stream << setw(6) << setprecision(2) << fixed << right << 1.0;
		stream << setw(12) << right << a.element();
		stream << setw(2) << "  ";
		stream << setw(5) << right << a.idatm_type_unmask();
		stream << setw(5) << right << a.gaff_type();
		stringstream ss;
		ss << " aps={";
		for (auto it = a.__aps.begin(); it != a.__aps.end(); ++it) {
			auto it_plus_one = it;
			ss << "{" << it->first << "," << it->second << "}" 
				<< (++it_plus_one == a.__aps.end() ? "" : ",");
		}
		stream << ss.str() << "}"; 
		ss.str("");
		ss << " prop={";
		for (auto it = a.__smiles_prop.begin(); it != a.__smiles_prop.end(); ++it) {
			auto it_plus_one = it;
			ss << "{" << it->first << "," << it->second << "}" 
				<< (++it_plus_one == a.__smiles_prop.end() ? "" : ",");
		}
		stream << ss.str() << "}" << endl; 
		return stream;
	}

	ostream& operator<< (ostream& stream, const Atom::Set &atoms) {
		for (auto &patom : atoms) stream << *patom;
		return stream;
	}
	ostream& operator<< (ostream& stream, const Atom::Vec &atoms) {
		for (auto &patom : atoms) stream << *patom;
		return stream;
	}

	BondSet get_bonds_in(const Atom::Set &atoms, bool in) {
		BondSet bonds;
		for (auto &patom1 : atoms) {
			for (auto &pbond : patom1->get_bonds()) {
				Atom &atom2 = pbond->second_atom(*patom1);
				if ((in || !atoms.count(&atom2)) 
					&& (!in || atoms.count(&atom2))) {
					bonds.insert(pbond);
				}
			}
		}
		return bonds;
	}

	BondSet get_bonds_in(const Atom::Vec &atoms, bool in) {
		Atom::Set a(atoms.begin(), atoms.end());
		return get_bonds_in(a, in);
	}

	Atom::Graph Atom::create_graph(const Atom::Vec &atoms) {
		return Atom::Graph(atoms, true);
	}

	Atom::Graph Atom::create_graph(const Atom::Set &atoms) {
		return Atom::Graph(atoms, true);
	}

	int Atom::get_num_hydrogens() const {
		const Atom &atom = *this;
		int num_h = 0;
		for (auto &atom2 : atom)
			if (atom2.element() == Element::H) 
				num_h++;
		return num_h;
	}
	Bond& Atom::connect(Atom &a2) {
		Atom &a1 = *this;
		if (!a1.is_adjacent(a2)) {
			dbgmsg("connecting atoms " << a1.atom_number() << " and " << a2.atom_number());
			a1.add(&a2);
			a2.add(&a1);
			a1.insert_bond(a2, new Bond(&a1, &a2)); // insert if not exists
			return *a2.insert_bond(a1, a1.get_shared_ptr_bond(a2)); // insert if not exists
		}
		return a1.get_bond(a2); // if already connected return existing bond
	}
	void Atom::set_members(const string &str) {
		// set atomic penalty scores
		boost::smatch m;
		dbgmsg("set members for atom " << this->atom_name() << " " << this->atom_number());
		if (boost::regex_search(str, m, boost::regex("aps=(\\S+)"))) {
			if (m[1].matched) {
				string aps_str = m[1].str();
				dbgmsg("aps_str = " << aps_str);
				boost::match_results<string::const_iterator> matches;
				string::const_iterator begin = aps_str.begin(), end = aps_str.end();
				while (boost::regex_search(begin, end, matches, boost::regex("\\{(\\d+,\\d+)\\}"))) {
					string s(matches[1].first, matches[1].second);
					auto vec = help::ssplit(s, ",");
					int val = stoi(vec[0]);
					int aps = stoi(vec[1]);
					__aps[val] = aps;
					begin = matches[1].second;
				}
			}
		}
		// set properties
		if (boost::regex_search(str, m, boost::regex("prop=(\\S+)"))) {
			if (m[1].matched) {
				string prop_str = m[1].str();
				dbgmsg("prop = " << prop_str);
				boost::match_results<string::const_iterator> matches;
				string::const_iterator begin = prop_str.begin(), end = prop_str.end();
				while (boost::regex_search(begin, end, matches, boost::regex("\\{(\\w+,\\d+)\\}"))) {
					string s(matches[1].first, matches[1].second);
					auto vec = help::ssplit(s, ",");
					string val = vec[0];
					int count = stoi(vec[1]);
					__smiles_prop[val] = count;
					begin = matches[1].second;
				}
			}
		}
		// set gaff type
		if (boost::regex_search(str, m, boost::regex("gaff=([^,$]+)"))) {
			if (m[1].matched) {
				dbgmsg(m[1].str());
				set_gaff_type(m[1].str());
			}
		}
		// set idatm type
		if (boost::regex_search(str, m, boost::regex("idatm=([^,$]+)"))) {
			if (m[1].matched) {
				dbgmsg(m[1].str());
				set_idatm_type(m[1].str());
			}
		}
	}

	bool Atom::compatible(const Atom &other) const {
		/* this is "smiles" atom and other is real atom */
		if (!boost::regex_search(other.get_label(), boost::regex(this->get_label())))
			return false;
		for (auto &kv : this->__smiles_prop) {
			auto &p = kv.first;
			auto &cnt = kv.second;
			int other_cnt = other.compute_num_property(p);
			if (!((cnt == -1 && other_cnt > 0) || cnt == other_cnt))
				return false;
		}
		return true;
	}
	int Atom::compute_num_property(const string &prop) const {
		if (prop == "B") return this->size();
		else if (prop == "H") return this->get_num_hydrogens();
		else if (prop == "sb" || prop == "db" || prop == "tb" ||
			//~ prop == "SB" || prop == "DB" || prop == "TB" || prop == "AB") 
			prop == "SB" || prop == "DB" || prop == "DL") 
			return this->get_num_bond_with_bond_gaff_type(prop);
		else if (prop.substr(0,2) == "AR" || prop.substr(0,2) == "RG"
			 || prop.substr(0,2) == "NG" || prop.substr(0,2) == "ag")
			return this->get_num_property(prop);
		else
			throw Error("Invalid property: " + prop);
	}
	int Atom::get_num_bond_with_bond_gaff_type(const string &prop) const {
		const Atom &atom = *this;
		int num_bwp = 0;
		for (auto &pbond : atom.get_bonds())
			if (pbond->get_bond_gaff_type() == prop) 
				num_bwp++;
		return num_bwp;
	}

        double molib::Atom::compute_rmsd(const Graph& g1, const Graph& g2) {
                Atom::Graph::Matches m = g1.match(g2);

                // try calculating rmsd of each mapping of molecule to molecule ...
                double min_sum_squared = HUGE_VAL;

                for (auto &mv : m) {
                        auto &vertices1 = mv.first;
                        auto &vertices2 = mv.second;
                        // match must stretch over the whole graph g1 (and g2)
                        if (vertices2.size() != g1.size()) 
                                throw Error("die : RMSD calculation aborted due to imperfect match"); 
                        double sum_squared = 0;
                        for (size_t i = 0; i < vertices2.size(); ++i) {
                                Atom &a1 = g1[vertices1[i]];
                                Atom &a2 = g2[vertices2[i]];
                                sum_squared += a1.crd().distance_sq(a2.crd());
                        }
                        if (sum_squared < min_sum_squared)
                                min_sum_squared = sum_squared;
                }

                return sqrt(min_sum_squared / g1.size());	
        }

        double compute_rmsd_vina_sq(const Atom::Vec &crds1, const Atom::Vec &crds2) {

                double distance_sum = 0.0;

                for (const auto &a : crds1) {
                        double min_dist = std::numeric_limits<double>::max();
                        for (const auto &b : crds2) {

                                if (a->idatm_type() != b->idatm_type()) {
                                        continue;
                                }

                                double temp = b->crd().distance_sq(a->crd());
                                min_dist = std::min(min_dist, temp);
                        }
                        distance_sum += min_dist;
                }

                return distance_sum / ( (crds1.size() + crds2.size()) / 2 );
        }

        double compute_rmsd_vina(const Atom::Vec &crds1, const Atom::Vec &crds2) {
                double rmsd_sq = compute_rmsd_vina_sq(crds1, crds2);
                return std::sqrt(rmsd_sq);
        }


};
}
