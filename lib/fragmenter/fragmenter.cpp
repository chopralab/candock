#include "fragmenter.hpp"
#include "unique.hpp"
#include <queue>
#include <iterator>
#include "pdbreader/molecule.hpp"
#include "helper/help.hpp"
#include "graph/graph.hpp"

namespace Molib {
	ostream& operator<<(ostream& os, const AtomMatch& atom_match) {
		for (auto &kv : atom_match) {
			os << "smiles atom number = " << kv.first 
				<< " matched atom = " << *kv.second << endl;
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

	//~ Fragmenter::Fragmenter(const AtomVec &atoms) : __mol_graph(create_graph(atoms)), 
		//~ __min_rigid_atoms(5), __min_gaff_group_atoms(3) {}
	Fragmenter::Fragmenter(const AtomVec &atoms) : __atoms(atoms), 
		__min_rigid_atoms(5), __min_gaff_group_atoms(3) {}
	//~ Fragmenter::Fragmenter(Model &model) : __mol_graph(create_graph(model)), 
		//~ __min_rigid_atoms(5), __min_gaff_group_atoms(3) {}
	//~ Fragmenter::Fragmenter(Residue &residue) : __mol_graph(create_graph(residue)), 
		//~ __min_rigid_atoms(5), __min_gaff_group_atoms(3) {}
	//~ Rings Fragmenter::__bond_rings_to_atom_rings(const set<BondSet> &brings) {
		//~ Rings arings;
		//~ for (auto &bring : brings) {
			//~ Ring aring;
			//~ for (auto &pbond : bring) {
				//~ aring.insert(&pbond->first_atom());
				//~ aring.insert(&pbond->second_atom());
			//~ }
			//~ arings.insert(aring);
		//~ }
		//~ return arings;
	//~ }
	
	//~ Rings Fragmenter::identify_fused_rings() {
		//~ set<BondSet> brings = __mol_graph.find_fused_rings();
		//~ return __bond_rings_to_atom_rings(brings);
		//~ // return __mol_graph.find_fused_rings();
	//~ }
	Rings Fragmenter::identify_fused_rings() {
		return create_graph(__atoms).find_fused_rings();
	}
	//~ Rings Fragmenter::identify_rings() {
		//~ set<BondSet> brings = __mol_graph.find_rings();
		//~ return __bond_rings_to_atom_rings(brings);
		//~ // return __mol_graph.find_rings();
	//~ }
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
	//~ void Fragmenter::__merge_small_fragments_with_rings(AtomSet &ring, int frag_size) {
		//~ for (auto &pa : ring) {
			//~ // extend from each atom of a ring by max frag_size
			//~ queue<Atom*> q;
			//~ AtomSet visited;
			//~ q.push(pa);
			//~ // visited.insert(pa);
			//~ while (!q.empty()) {
				//~ Atom &atom = *q.front(); q.pop();
				//~ visited.insert(&atom);
				//~ // for (auto &adj_a : *t) {
				//~ for (auto &bond : atom) {
					//~ // if (!ring.count(&adj_a) && !visited.count(&adj_a)) {
					//~ if (!ring.count(&bond.second_atom()) 
						//~ && !visited.count(&bond.second_atom())) {
						//~ q.push(&bond.second_atom());
						//~ // visited.insert(&bond.second_atom());
					//~ }
				//~ }
				//~ if (visited.size() > frag_size) goto skip_insert;
			//~ }
			//~ ring.insert(visited.begin(), visited.end());
			//~ skip_insert:
			//~ ;
		//~ }
	//~ }
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
				//~ for (auto &bond : atom) {
					if (!ring.count(&adj_a) && !visited.count(&adj_a)) {
					//~ if (!ring.count(&bond.second_atom()) 
						//~ && !visited.count(&bond.second_atom())) {
						//~ q.push(&bond.second_atom());
						q.push(&adj_a);
						//~ visited.insert(&bond.second_atom());
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
					//~ for (auto &bond : current) {
						//~ auto &bondee = bond.second_atom();
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
							//~ q.push(&bond.second_atom());
							q.push(&bondee);
						}
					}
				}
			}
		}
	}
	//~ Bonds Fragmenter::identify_rotatable(const help::rotatable_bonds &bonds) {
		//~ Bonds rotatable;
		//~ for (auto &kv : bonds) {
			//~ string bond_type = kv.first;
			//~ help::bdef bond_data = kv.second;
			//~ AtomTags atom_tags = create_atom_tags(get<0>(bond_data));
			//~ Glib::Graph<AtomTag> g(atom_tags, true);
			//~ dbgmsg(g);
			//~ Glib::Graph<AtomTag>::Matches matches = g.match(__mol_graph);
			//~ for (auto &match : matches) { // apply fragment name + id to all found fragments
				//~ AtomVec atoms;
				//~ for(auto &nid : match.second) atoms.push_back(&__mol_graph[nid]);
				//~ Atom &a1 = (atoms[0]->atom_number() < atoms[1]->atom_number() ? *atoms[0] : *atoms[1]);
				//~ Atom &a2 = (atoms[0]->atom_number() < atoms[1]->atom_number() ? *atoms[1] : *atoms[0]);
				//~ rotatable[BondKey(&a1, &a2)] = Bond(&a1, &a2, bond_type);  // overwrite existing bond type with more specialized one
			//~ }
//~ #ifndef NDEBUG
			//~ dbgmsg("matches for " << bond_type);
			//~ for (auto &match : matches) {
				//~ dbgmsg("match = ");
				//~ int i=0;
				//~ for (auto &node_id1 : match.first) {
					//~ Glib::node_id node_id2 = match.second[i++];
					//~ dbgmsg("[" << g[node_id1].get_label() << "," << __mol_graph[node_id2].get_label());
				//~ }
			//~ }
			//~ dbgmsg("after matches");
//~ #endif
		//~ }
		//~ // rotatable bond inside a ring is changed back to non-rotatable
		//~ Rings rings = identify_fused_rings();
		//~ for (auto &ring : rings) {
			//~ for (auto &pa1 : ring) {
				//~ for (auto &pa2 : ring) {
					//~ rotatable.erase(BondKey(pa1, pa2));
				//~ }
			//~ }
		//~ }
//~ #ifndef NDEBUG
		//~ dbgmsg("rotatable bonds:");
		//~ for (auto &kv : rotatable) {
			//~ Bond &b = kv.second;
			//~ dbgmsg(b.first_atom().idatm_type() << "_" << b.first_atom().atom_number() 
				//~ << "-" << b.second_atom().idatm_type() << "_" << b.second_atom().atom_number() 
				//~ << " " << b.get_type());
		//~ }
//~ #endif
		//~ dbgmsg("exiting identify_rotatable");
		//~ return rotatable;
	//~ }
	
	//~ void Fragmenter::apply_rule(const AtomMatch &patt, const vector<string> &bond_rules) {
		//~ for (auto &bond_rule : bond_rules) {
			//~ const vector<string> s = help::ssplit(bond_rule, ":", true);
			//~ if (s.size() < 3) throw Error("die : wrong bond rule");
			//~ int atom_number1 = stoi(s[0]);
			//~ int atom_number2 = stoi(s[1]);
			//~ string bond_type = s[2];
			//~ Atom &a1 = patt.at(atom_number1);
			//~ Atom &a2 = patt.at(atom_number2);
			//~ Bond &bond12 = a1[&a2];
			//~ Bond &bond21 = a2[&a1];
			//~ if (bond12.get_rotatable().empty()) { // don't overwrite
				//~ bond12.set_rotatable(bond_type);
				//~ bond21.set_rotatable(bond_type);
			//~ }
		//~ }
	//~ }
	//~ void Fragmenter::apply_rule(const AtomMatch &m, const vector<string> &rules,
		//~ AtomSet &visited_atoms, BondSet &visited_bonds) {
		//~ for (auto &rule : rules) {
			//~ auto vec1 = help::ssplit(rule, ":", true);
			//~ auto vec2 = help::ssplit(vec1[0], ",", true);
			//~ if (vec2.size() == 1) { // atom rule
				//~ int atom_number = stoi(vec2[0]);
				//~ Atom &atom = m.at(atom_number);
				//~ if (!visited_atoms.count(&atom)) {
					//~ atom.set_members(vec1[1]);
					//~ visited_atoms.insert(&atom);
				//~ }
			//~ } else if (vec2.size() == 2) { // bond rule
				//~ int atom_number1 = stoi(vec2[0]);
				//~ int atom_number2 = stoi(vec2[1]);
				//~ Atom &atom1 = m.at(atom_number1);
				//~ Atom &atom2 = m.at(atom_number2);
				//~ Bond &bond12 = atom1[&atom2];
				//~ Bond &bond21 = atom2[&atom1];
				//~ if (!visited_bonds.count(&bond12)) {
					//~ bond12.set_members(vec1[1]);
					//~ bond21.set_members(vec1[1]);
					//~ visited_bonds.insert(&bond12);
					//~ visited_bonds.insert(&bond21);
				//~ }
			//~ }
		//~ }
	//~ }
	//~ 
	void Fragmenter::apply_rule(const AtomMatch &m, 
		const vector<string> &rules, AtomSet &visited) { // atom rule
		for (auto &rule : rules) {
			dbgmsg(rule);
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
			//~ Bond &bond12 = atom1[&atom2];
			Bond &bond12 = atom1.get_bond(atom2);
			//~ Bond &bond21 = atom2[&atom1];
			if (!visited.count(&bond12)) {
				bond12.set_members(vec1[1]);
				//~ bond21.set_members(vec1[1]);
				visited.insert(&bond12);
				//~ visited.insert(&bond21);
			}
		}
	}
	void Fragmenter::substitute_bonds(const help::rename_rules &rrules) {
		BondSet visited;
		Fragmenter &fragmenter = *this;
		for (auto &p : rrules) {
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
		}
	}
	//~ Fragmenter::AtomMatchVec Fragmenter::grep(const help::smiles &smi) {
		//~ AtomMatchVec pv;
		//~ MolGraph g = create_graph(smi);
		//~ MolGraph::Matches matches = g.match(__mol_graph);
		//~ for (auto &match : matches) {
			//~ Pattern patt;
			//~ for(int i = 0; i < match.first.size(); i++) {
				//~ auto &nid1 = match.first[i];
				//~ auto &nid2 = match.second[i];
				//~ Bond &bond1 = g[nid1];
				//~ Bond &bond2 = __mol_graph[nid2];
				//~ Atom &a11 = bond1.first_atom();
				//~ Atom &a12 = bond1.second_atom();
				//~ Atom &a21 = bond2.first_atom();
				//~ Atom &a22 = bond2.second_atom();
				//~ patt[a11.atom_number()] = &a21;
				//~ patt[a12.atom_number()] = &a22;
			//~ }
			//~ pv.push_back(patt);
		//~ }
		//~ return p;
	//~ }
	//~ AtomMatchVec Fragmenter::grep(const help::smiles &smi) {
		//~ AtomMatchVec mvec;
		//~ MolGraph smiles_graph = create_graph(smi);
		//~ MolGraph::Matches matches = smiles_graph.match(__mol_graph);
		//~ for (auto &match : matches) {
			//~ AtomMatch m;
			//~ for(int i = 0; i < match.first.size(); ++i) {
				//~ auto &nid1 = match.first[i];
				//~ auto &nid2 = match.second[i];
				//~ Bond &bond1 = smiles_graph[nid1];
				//~ Bond &bond2 = __mol_graph[nid2];
				//~ Atom &a11 = bond1.first_atom();
				//~ Atom &a12 = bond1.second_atom();
				//~ Atom &a21 = bond2.first_atom();
				//~ Atom &a22 = bond2.second_atom();
				//~ m[a11.atom_number()] = &a21;
				//~ m[a12.atom_number()] = &a22;
			//~ }
			//~ mvec.push_back(m);
		//~ }
		//~ return mvec;
	//~ }
	//~ AtomMatchVec Fragmenter::grep(const help::smiles &smi) {
		//~ AtomMatchVec mvec;
		//~ BondGraph smiles_graph = create_graph(smi);
		//~ BondGraph bond_graph = create_graph(get_bonds_in(__atoms));
		//~ BondGraph::Matches matches = smiles_graph.match(bond_graph);
		//~ for (auto &match : matches) {
			//~ AtomMatch m;
			//~ for(int i = 0; i < match.first.size(); ++i) {
				//~ auto &nid1 = match.first[i];
				//~ auto &nid2 = match.second[i];
				//~ Bond &bond1 = smiles_graph[nid1];
				//~ Bond &bond2 = bond_graph[nid2];
				//~ if (!m.count(bond1.atom1().atom_number()) 
					//~ && !m.count(bond1.atom2().atom_number())) {
				//~ }
					//~ 
				//~ m[bond1.atom1().atom_number()] = &bond2.atom1();
				//~ m[bond1.atom2().atom_number()] = &bond2.atom2();
			//~ }
			//~ mvec.push_back(m);
		//~ }
		//~ return mvec;
	//~ }
	AtomMatchVec Fragmenter::grep(const help::smiles &smi) {
		AtomMatchVec mvec;
		BondGraph smiles_graph = create_graph(smi);
		BondGraph bond_graph = create_graph(get_bonds_in(__atoms));
		BondGraph::Matches matches = smiles_graph.match(bond_graph);
		for (auto &match : matches) {
			map<Bond*, Bond*> bond_match;
			for (int i = 0; i < match.first.size(); ++i) {
				auto &nid1 = match.first[i];
				auto &nid2 = match.second[i];
				Bond &bond1 = smiles_graph[nid1];
				Bond &bond2 = bond_graph[nid2];
				bond_match[&bond1] = &bond2;
			}
			dbgmsg("in grep : bond match is " << endl << bond_match);
			try {
				AtomMatch m = __convert_to_atom_match(bond_match);
				mvec.push_back(m);
				dbgmsg("in grep : atom match (try 1) is " << endl << m);
			} catch(Error &e) {
				try {
					AtomMatch m = __convert_to_atom_match(bond_match, true);
					mvec.push_back(m);
					dbgmsg("in grep : atom match (try 2) is " << endl << m);
				} catch(Error &e) {}
			}
		}
		return mvec;
	}
	//~ Fragmenter::AtomMatch Fragmenter::__convert_to_atom_match(const map<Bond*, Bond*> &bond_match) {
		//~ AtomMatch amatch;
		//~ map<Bond*, Bond*> bmatch = bond_match;
		//~ bool repeat = false;
		//~ while (!bmatch.empty()) {
			//~ for (auto &kv : bmatch) {
				//~ Bond &bond1 = *kv.first;
				//~ Bond &bond2 = *kv.second;
				//~ Atom &atom11 = bond1.atom1();
				//~ Atom &atom12 = bond1.atom2();
				//~ Atom &atom21 = bond2.atom1();
				//~ Atom &atom22 = bond2.atom2();
				//~ if (!repeat) {
				//~ if (atom11.compatible(atom21)) {
					//~ amatch[atom11.atom_number()] = &atom21;
					//~ amatch[atom12.atom_number()] = &atom22;
				//~ }
				//~ bmatch.erase(&bond1);
				//~ break;
			//~ }
			//~ while (!bmatch.empty()) {
				//~ int sz = bmatch.size();
				//~ for (auto &kv : bmatch) {
					//~ Bond &bond1 = *kv.first;
					//~ Bond &bond2 = *kv.second;
					//~ int anum11 = bond1.atom1().atom_number();
					//~ int anum12 = bond1.atom2().atom_number();
					//~ Atom &atom21 = bond2.atom1();
					//~ Atom &atom22 = bond2.atom2();
					//~ if (amatch1.count(anum11)) {
						//~ amatch1[anum12] = (amatch1.at(anum11) == &atom21 ? &atom22 : &atom21);
						//~ bmatch.erase(&bond1);
					//~ } else if (amatch1.count(&anum12)) {
						//~ amatch1[anum11] = (amatch1.at(anum12) == &atom21 ? &atom22 : &atom21);
						//~ bmatch.erase(&bond1)
					//~ }
				//~ }
				//~ if (bmatch.size() == sz) { 
					//~ bmatch = bond_match; 
					//~ amatch.clear();
					//~ repeat = true; 
					//~ break; 
				//~ }
			//~ }
		//~ }
		//~ return amatch;
	//~ }
	//~ Fragmenter::AtomMatch Fragmenter::__convert_to_atom_match(const map<const Bond*, Bond*> &bond_match) {
		//~ AtomMatch amatch;
		//~ queue<Bond*> q;
		//~ q.push(bond_match.begin()->first);
		//~ while (!q.empty()) {
			//~ Bond &bond1 = *q.front(); q.pop();		
			//~ Bond &bond2 = bond_match.at(&bond1);
			//~ int anum11 = bond1.atom1().atom_number();
			//~ int anum12 = bond1.atom2().atom_number();
			//~ Atom &atom21 = bond2.atom1();
			//~ Atom &atom22 = bond2.atom2();
//~ 
			//~ if (amatch[anum11] == &atom21) amatch[anum12] = &atom22;
			//~ if (amatch[anum11] == &atom22) amatch[anum12] = &atom21;
			//~ if (amatch[anum12] == &atom21) amatch[anum11] = &atom22;
			//~ if (amatch[anum12] == &atom22) amatch[anum11] = &atom21;
			//~ visited.insert(&bond1);
			//~ for (auto &pbond : bond1)
				//~ if (!visited.count(&pbond) && bond_match.count(pbond))
					//~ q.push(pbond);
		//~ }
		//~ return amatch;
	//~ }
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
				if (amatch[anum11] == &atom21) amatch[anum12] = &atom22;
				else if (amatch[anum11] == &atom22) amatch[anum12] = &atom21;
				else if (amatch[anum12] == &atom21) amatch[anum11] = &atom22;
				else if (amatch[anum12] == &atom22) amatch[anum11] = &atom21;
				else throw Error("exception : cannot convert from bond match to atom match");
			}
			visited.insert(&bond1);
			for (auto &bond : bond1)
				if (!visited.count(&bond) && bond_match.count(&bond))
					q.push(&bond);
		}
		return amatch;
	}
	//~ 
	//~ Fragmenter::AtomToRename Fragmenter::rename(const help::rename_rules &r) {
		//~ dbgmsg("entering rename");
		//~ AtomToRename rename_atom;
		//~ for (auto &kv : r) {
			//~ help::smiles v_smiles = kv.first;
			//~ help::rname m_rn = kv.second;
			//~ dbgmsg("before create graph");
			//~ AtomTags atom_tags = create_atom_tags(v_smiles);
			//~ Glib::Graph<AtomTag> g(atom_tags, true);
			//~ dbgmsg("graph is " << g);
			//~ MolGraph::Matches matches = g.match(__mol_graph);
			//~ for (auto &match : matches) {
				//~ for(int i = 0; i < match.first.size(); i++) {
					//~ auto &nid1 = match.first[i];
					//~ auto &nid2 = match.second[i];
					//~ if (m_rn.count(nid1 + 1)) {
						//~ rename_atom[&__mol_graph[nid2]] = m_rn.at(nid1 + 1); // overwrite with more specific type
					//~ }
				//~ }
			//~ }
		//~ }
		//~ return rename_atom;
	//~ }
	
	//~ Fragments Fragmenter::identify_overlapping_rigid_segments(const AtomSet &atoms) {
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
					//~ for (auto &bond : current)
						//~ if (!visited.count(&bond.second_atom()) 
							//~ && !bond.is_rotatable())
							//~ q.push(&bond.second_atom());
					for (auto &pbond : current.get_bonds())
						if (!visited.count(&pbond->second_atom(current)) 
							&& !pbond->is_rotatable())
							q.push(&pbond->second_atom(current));
				}
				// insert atoms that are connected to the segment 
				// by one rotatable bond
				//~ for (auto &pbond : get_bonds_in(seg, false)) {
					//~ Bond &bond = *pbond;
					//~ if (bond.is_rotatable())
						//~ seg.insert(&bond.second_atom());
				//~ }
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

