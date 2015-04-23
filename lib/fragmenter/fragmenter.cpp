#include "fragmenter.hpp"
#include "unique.hpp"
#include <queue>
#include <iterator>
#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"
#include "graph/graph.hpp"

namespace Molib {
	ostream& operator<<(ostream& os, const AtomMatchVec& atom_match_vec) {
		for (auto &atom_match : atom_match_vec) {
			os << atom_match << endl;
		}
		return os;
	}

	ostream& operator<<(ostream& os, const AtomMatch& atom_match) {
		for (auto &kv : atom_match) {
			os << "[" << kv.first 
				<< "-->" << kv.second->atom_name() << "] ";
		}
		return os;
	}

	ostream& operator<<(ostream& os, const map<Bond*, Bond*>& bond_match) {
		for (auto &kv : bond_match) {
			os << "smiles bond = " << *kv.first 
				<< " matched molecule bond = " << *kv.second << endl;
		}
		return os;
	}

	ostream& operator<<(ostream& os, const Fragments& fragments)	{
		for (auto &fragment : fragments) {
			os << "fragment " << fragment.first << " : " << endl;
			for (auto &atoms : fragment.second) {
				for (auto &pa : atoms) {
					os << "        " << *pa;
				}
				os << "-------------------";
			}
		}
		return os;
	}

	Fragmenter::Fragmenter(const AtomVec &atoms) : __atoms(atoms), 
		__min_rigid_atoms(5), __min_gaff_group_atoms(3) {}

	Rings Fragmenter::identify_fused_rings() {
		return create_graph(__atoms).find_fused_rings();
	}

	Rings Fragmenter::identify_rings() {
		return create_graph(__atoms).find_rings();
	}

	Fragments Fragmenter::identify_seeds(const Fragments &rigid, Unique &u) {
		/* Seeds are essentially rigid segments with > 4 atoms (e.g., rings, 
		 * non-ring rigid fragments, etc.) to which 
		 * small groups with < 4 atoms were added (e.g., -NO2, -CH2-CH3, etc.).
		 * Note that there can be a rotatable bond between rigid segment and 
		 * a small group. This introduces a small error in the maximum clique
		 * part which is rigid only.
		 * 
		 */
		Fragments seeds;
		for (auto &kv : rigid) {
			string nm = kv.first;
			const set<AtomSet> &rigids_with_same_name = kv.second;
			if (rigids_with_same_name.empty() || rigids_with_same_name.size() > 1)
				throw Error("die : there should be only one rigid segment per one name");
			const AtomSet &rigid_segment = *rigids_with_same_name.begin();
			if (rigid_segment.size() > __min_rigid_atoms) {
				AtomSet a(rigid_segment.begin(), rigid_segment.end());
				__merge_small_fragments_with_rings(a, __min_gaff_group_atoms);
				seeds[help::to_string(u.get_seed_id(a))].insert(a);
			}
		}
		dbgmsg("------------- SEEDS -------------------" << endl << seeds);
		return seeds;
	}

	void Fragmenter::__merge_small_fragments_with_rings(AtomSet &ring, int frag_size) {
		for (auto &pa : ring) {
			// extend from each atom of a ring by max frag_size
			queue<Atom*> q;
			AtomSet visited;
			q.push(pa);
			while (!q.empty()) {
				Atom &atom = *q.front(); q.pop();
				visited.insert(&atom);
				for (auto &adj_a : atom) {
					if (!ring.count(&adj_a) && !visited.count(&adj_a)) {
						q.push(&adj_a);
					}
				}
				if (visited.size() > frag_size) goto skip_insert;
			}
			ring.insert(visited.begin(), visited.end());
			skip_insert:
			;
		}
	}
	void Fragmenter::flip_conjugated_gaff_types(const AtomVec &atoms) {
		AtomSet visited;
		for (auto &patom : atoms) {
			if (!visited.count(patom) && (help::gaff_group_1.count(patom->gaff_type()) 
				|| help::gaff_group_2.count(patom->gaff_type()))) {
				// visit bonds in every direction from current atom
				// and flip conjugated gaff types
				queue<Atom*> q;
				q.push(patom);
				while (!q.empty()) {
					Atom &current = *q.front(); q.pop();
					visited.insert(&current);
					for (auto &bondee : current) {
						Bond &bond = current.get_bond(bondee);
						if (!visited.count(&bondee) && (help::gaff_group_1.count(bondee.gaff_type()) 
							|| help::gaff_group_2.count(bondee.gaff_type()))) {
							if (bond.is_double() || bond.is_triple() || bond.is_aromatic()) {
								if (help::gaff_group_1.count(current.gaff_type()) 
									&& help::gaff_group_1.count(bondee.gaff_type())
									|| help::gaff_group_2.count(current.gaff_type()) 
									&& help::gaff_group_2.count(bondee.gaff_type())) {
									dbgmsg("bondee gaff type (case 1) = "
										<< bondee.gaff_type());
									bondee.set_gaff_type(help::gaff_flip.at(bondee.gaff_type()));
									visited.insert(&bondee);
								}
							} else if (bond.is_single()) {
								if (help::gaff_group_1.count(current.gaff_type()) 
									&& help::gaff_group_2.count(bondee.gaff_type())
									|| help::gaff_group_2.count(current.gaff_type()) 
									&& help::gaff_group_1.count(bondee.gaff_type())) {
									dbgmsg("bondee gaff type (case 2) = "
										<< bondee.gaff_type());
									bondee.set_gaff_type(help::gaff_flip.at(bondee.gaff_type()));
									visited.insert(&bondee);
								}
							}
							q.push(&bondee);
						}
					}
				}
			}
		}
	}

	void Fragmenter::apply_rule(const AtomMatch &m, 
		const vector<string> &rules, AtomSet &visited) { // atom rule
		for (auto &rule : rules) {
			dbgmsg("atom_match = " << m << " rule = " << rule);
			auto vec1 = help::ssplit(rule, ":", true);
			assert(vec1.size() == 2);
			dbgmsg(vec1[0] << " " << vec1[1]);
			int atom_number = stoi(vec1[0]);
			Atom &atom = *m.at(atom_number);
			if (!visited.count(&atom)) {
				atom.set_members(vec1[1]);
				visited.insert(&atom);
			}
		}
	}
	
	void Fragmenter::apply_rule(const AtomMatch &m, 
		const vector<string> &rules, BondSet &visited) { // bond rule
		for (auto &rule : rules) {
			auto vec1 = help::ssplit(rule, ":", true);
			auto vec2 = help::ssplit(vec1[0], ",", true);
			assert(vec1.size() == 2);
			assert(vec2.size() == 2);
			dbgmsg(vec1[0] << " " << vec1[1]);
			dbgmsg(vec2[0] << " " << vec2[1]);
			int atom_number1 = stoi(vec2[0]);
			int atom_number2 = stoi(vec2[1]);
			Atom &atom1 = *m.at(atom_number1);
			Atom &atom2 = *m.at(atom_number2);
			Bond &bond12 = atom1.get_bond(atom2);
			if (!visited.count(&bond12)) {
				bond12.set_members(vec1[1]);
				visited.insert(&bond12);
			}
		}
	}

	void Fragmenter::substitute_bonds(const help::rename_rules &rrules) {
		dbgmsg("starting substitute_bonds");
		BondSet visited;
		Fragmenter &fragmenter = *this;
		for (auto &p : rrules) {
			dbgmsg("grepping for rule = " << p);
			AtomMatchVec atom_matches = fragmenter.grep(p.pattern);
			for (auto &m : atom_matches)
				fragmenter.apply_rule(m, p.rule, visited);
		}
	}

	void Fragmenter::substitute_atoms(const help::rename_rules &rrules) {
		AtomSet visited;
		Fragmenter &fragmenter = *this;
		for (auto &p : rrules) {
			AtomMatchVec atom_matches = fragmenter.grep(p.pattern);
			for (auto &m : atom_matches)
				fragmenter.apply_rule(m, p.rule, visited);
			dbgmsg("RENAME RULE : " << p << endl
				<< "MATCHES : " << atom_matches);

		}
	}

	AtomMatchVec Fragmenter::grep(const help::smiles &smi) {
		AtomMatchVec mvec;
		BondGraph smiles_graph = create_graph(smi);
		BondGraph bond_graph = create_graph(get_bonds_in(__atoms));
		BondGraph::Matches matches = smiles_graph.match(bond_graph);
		dbgmsg("NUM MATCHES FOUND = " << matches.size());
		for (auto &match : matches) {
			map<Bond*, Bond*> bond_match;
			dbgmsg("MATCH SIZE = " << match.first.size());
			for (int i = 0; i < match.first.size(); ++i) {
				auto &nid1 = match.first[i];
				auto &nid2 = match.second[i];
				Bond &bond1 = smiles_graph[nid1];
				Bond &bond2 = bond_graph[nid2];
				bond_match[&bond1] = &bond2;
			}
			dbgmsg("in grep : bond match is " << endl << bond_match);
			for (bool reverse : {false, true}) { // two tries : forward and reverse
				try {
					AtomMatch m = __convert_to_atom_match(bond_match, reverse);
					mvec.push_back(m);
					dbgmsg("in grep : atom match (try " << reverse << ") is " << endl << m);
				} catch(Error &e) {}
			}
		}
		return mvec;
	}

	AtomMatch Fragmenter::__convert_to_atom_match(
		const map<Bond*, Bond*> &bond_match, bool reverse) {
		AtomMatch amatch;
		BondSet visited;
		queue<Bond*> q;
		q.push(bond_match.begin()->first);
		while (!q.empty()) {
			Bond &bond1 = *q.front(); q.pop();		
			Atom &atom11 = bond1.atom1();
			Atom &atom12 = bond1.atom2();
			int anum11 = atom11.atom_number();
			int anum12 = atom12.atom_number();

			const Bond &bond2 = *bond_match.at(&bond1);
			Atom &atom21 = bond2.atom1();
			Atom &atom22 = bond2.atom2();

			if (visited.empty()) { // first bond-match considered
				if (!reverse) {
					if (atom11.compatible(atom21)) {
						amatch[anum11] = &atom21;
						amatch[anum12] = &atom22;
					} else throw Error("exception : trying with reversed first bond");
				} else {
					if (atom11.compatible(atom22)) {
						amatch[anum11] = &atom22;
						amatch[anum12] = &atom21;
					} else throw Error("exception : ran out of options");
				}
			} else { // not first bond
				if (amatch[anum11] == &atom21 && atom12.compatible(atom22)) amatch[anum12] = &atom22;
				else if (amatch[anum11] == &atom22 && atom12.compatible(atom21)) amatch[anum12] = &atom21;
				else if (amatch[anum12] == &atom21 && atom11.compatible(atom22)) amatch[anum11] = &atom22;
				else if (amatch[anum12] == &atom22 && atom11.compatible(atom21)) amatch[anum11] = &atom21;
				else throw Error("exception : cannot convert from bond match to atom match");
			}
			visited.insert(&bond1);
			for (auto &bond : bond1)
				if (!visited.count(&bond) && bond_match.count(&bond))
					q.push(&bond);
		}
		return amatch;
	}

	Fragments Fragmenter::identify_overlapping_rigid_segments(const AtomVec &atoms) {
		Fragments f;
		AtomSet visited;
		int i=0;
		for (auto &pa : atoms) {
			if (!visited.count(pa)) {
				// insert all atoms that are connected to this atom
				// by non-rotatable bond(s)
				AtomSet seg;
				queue<Atom*> q;
				q.push(pa);
				while (!q.empty()) {
					Atom &current = *q.front(); q.pop();
					seg.insert(&current);
					visited.insert(&current);
					for (auto &pbond : current.get_bonds())
						if (!visited.count(&pbond->second_atom(current)) 
							&& !pbond->is_rotatable())
							q.push(&pbond->second_atom(current));
				}
				// insert atoms that are connected to the segment 
				// by one rotatable bond
				for (auto &pbond : get_bonds_in(seg, false)) {
					Bond &bond = *pbond;
					if (bond.is_rotatable()) {
						seg.insert(&bond.atom1());
						seg.insert(&bond.atom2());
					}
				}
				// ensure segments have >= 3 atoms
				if (seg.size() >=3) {
					f[help::to_string(i++)].insert(seg);
				}
			}
		}
		dbgmsg("------------- RIGID SEGMENTS -------------------" << endl << f);
		return f;
	}
};

