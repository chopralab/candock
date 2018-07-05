#include "candock/linker/segment.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/help.hpp"

namespace candock {

namespace linker {
	ostream& operator<< (ostream& stream, const Segment& s) {
		stream << "Segment(" << s.__seed_id << ") : atom numbers = ";
		for (auto &pa : s.__atoms) stream << pa->atom_number() << " ";
		return stream;
	}

	Segment::Segment(const molib::Atom::Vec atoms, const int &seed_id, const Segment::Id idx) 
		: __atoms(atoms), __seed_id(seed_id), __id(idx), __join_atom(atoms.size(), false), 
		__common_atom(atoms.size(), false) { 

		for (size_t i = 0; i < atoms.size(); ++i) 
			__amap[atoms[i]] = i; 
	}

	//~ const molib::Atom& Segment::adjacent_in_segment(const molib::Atom &atom, 
		//~ const molib::Atom &forbidden) const { 
		//~ for (auto &adj : atom) {
			//~ if (&adj != &forbidden && has_atom(adj)) 
				//~ return adj; 
		//~ }
		//~ throw Error("die : couldn't find adjacent in segment");
	//~ }
//~ 
	int Segment::adjacent_in_segment(const molib::Atom &atom, 
		const molib::Atom &forbidden) const { 
		for (auto &adj : atom) {
			if (&adj != &forbidden && has_atom(adj)) 
				return get_idx(adj); 
		}
		throw Error("die : couldn't find adjacent in segment");
	}

	void Segment::set_bond(const Segment &other, molib::Atom &a1, molib::Atom &a2) { 
		__bond.insert({&other, molib::Bond(&a1, &a2, get_idx(a1), get_idx(a2))}); 
	}
	
	Segment::Graph Segment::create_graph(const molib::Molecule &molecule) {
		dbgmsg("Create segment graph ...");
		const molib::Model &model = molecule.first().first();
		vector<unique_ptr<Segment>> vertices;
		dbgmsg(model.get_rigid());
		int idx = 0;
		for (auto &fragment : model.get_rigid()) { // make vertices (segments) of a graph
			dbgmsg(fragment.get_all());
			auto all = fragment.get_all();
			molib::Atom::Vec fragatoms(all.begin(), all.end());
			//~ vertices.push_back(unique_ptr<Segment>(new Segment(fragatoms, fragment.get_seed_id())));
			vertices.push_back(unique_ptr<Segment>(new Segment(fragatoms, fragment.get_seed_id(), idx++)));
		}
		// connect segments
		for (size_t i = 0; i < vertices.size(); ++i) {
			Segment &s1 = *vertices[i];
			for (size_t j = i + 1; j < vertices.size(); ++j) {
				Segment &s2 = *vertices[j];
				auto inter = graph::intersection(s1.get_atoms(), s2.get_atoms());
				dbgmsg(s1.get_atoms().size() << " " << s2.get_atoms().size() << " " << inter.size());
				if (inter.size() == 1) {

					auto &atom = **inter.begin();
					s1.set_common_atom(atom);
					s2.set_common_atom(atom);
					dbgmsg("intersection size is one for segments " << s1 << " and " << s2);

				} else if (inter.size() == 2) {
					s1.add(&s2);
					s2.add(&s1);
					auto &atom1 = **inter.begin();
					auto &atom2 = **inter.rbegin();

					s1.set_common_atom(atom1);
					s2.set_common_atom(atom1);
					s1.set_common_atom(atom2);
					s2.set_common_atom(atom2);

					// determine which atom of bond is in s1 and which in s2
					int num_bonds = 0; 
					for (auto &adj : atom1) {
						if (s1.has_atom(adj))
							num_bonds++;
					}
					if (num_bonds == 1) {
						dbgmsg("atom " << atom1 << " belongs to segment " << s2);
						dbgmsg("atom " << atom2 << " belongs to segment " << s1);
						s1.set_bond(s2, atom2, atom1);
						s2.set_bond(s1, atom1, atom2);
						
						s1.set_join_atom(atom2);
						s2.set_join_atom(atom1);
					} else {
						dbgmsg("atom " << atom1 << " belongs to segment " << s1);
						dbgmsg("atom " << atom2 << " belongs to segment " << s2);
						s1.set_bond(s2, atom1, atom2);
						s2.set_bond(s1, atom2, atom1);

						s1.set_join_atom(atom1);
						s2.set_join_atom(atom2);

					}
				}
			}
		}

		const Segment::Paths paths = __find_paths(vertices);
		__init_max_linker_length(paths);
		__set_branching_rules(paths);

		return Segment::Graph(std::move(vertices), true, false);
	}

	Segment::Paths Segment::__find_paths(const vector<unique_ptr<Segment>> &segments) {
		/*
		 * Find ALL paths between ALL seed segments (even non-adjacent) 
		 * and seeds and leafs
		 */
		Segment::Paths paths;
		dbgmsg("find all paths in a graph");
		Segment::Set seeds, seeds_and_leafs;
		for (auto &pseg : segments) {
			auto &seg = *pseg;
			if (seg.is_seed())
				seeds.insert(&seg);
			if (seg.is_seed() || seg.is_leaf())
				seeds_and_leafs.insert(&seg);
		}
		set<Segment::ConstPair> visited;
		for (auto &pseg1 : seeds) {
			for (auto &pseg2 : seeds_and_leafs) {
				if (pseg1 != pseg2 && !visited.count({pseg1, pseg2})) {
					dbgmsg("finding path between " << *pseg1 << " and "
						<< *pseg2);
					visited.insert({pseg1, pseg2});
					visited.insert({pseg2, pseg1});
					Segment::Graph::Path path = graph::find_path(*pseg1, *pseg2);
					paths.insert({{pseg1, pseg2}, path });
					dbgmsg("path between first segment " << *pseg1
						<< " and second segment " << *pseg2);
				}
			}
		}
		return paths;
	}

	bool Segment::__link_adjacent(const Segment::Graph::Path &path) {
		/* Returns true if path has no seed segments along the way
		 * 
		 */
#ifndef NDEBUG
		for (auto it = path.begin(); it != path.end(); ++it) {
			dbgmsg(**it << " is seed = " << boolalpha << (*it)->is_seed());
		}
#endif
		for (auto it = path.begin() + 1; it != path.end() - 1; ++it) {
			dbgmsg(**it);
			if ((*it)->is_seed())
				return false;
		}
		return true;
	}

	void Segment::__init_max_linker_length(const Segment::Paths &paths) {
		for (auto &kv : paths) {
			Segment::Graph::Path path(kv.second.begin(), kv.second.end());
			__compute_max_linker_length(path);
		}
	}

	void Segment::__compute_max_linker_length(Segment::Graph::Path &path) {
		for (size_t j = 0; j < path.size() - 1; j++) {

			double d = 0.0;
			size_t i = j;
			
			molib::Atom *front_atom = &path[i]->get_bond(*path[i + 1]).atom1();
			
			for (; i < path.size() - 2; i += 2) {
				const molib::Bond &b1 = path[i]->get_bond(*path[i + 1]);
				const molib::Bond &b2 = path[i + 1]->get_bond(*path[i + 2]);
				path[j]->set_max_linker_length(*path[i + 1], d + b1.length());
				path[i + 1]->set_max_linker_length(*path[j], d + b1.length());
				
				d += front_atom->crd().distance(b2.atom2().crd());

				front_atom = &b2.atom2();
				
				path[j]->set_max_linker_length(*path[i + 2], d);
				path[i + 2]->set_max_linker_length(*path[j], d);

				dbgmsg("max_linker_length between " << *path[j] << " and " 
					<< *path[i + 1] << " = " << path[j]->get_max_linker_length(*path[i + 1]));
				dbgmsg("max_linker_length between " << *path[j] << " and " 
					<< *path[i + 2] << " = " << path[j]->get_max_linker_length(*path[i + 2]));
			}
			if (i < path.size() - 1) {
				const molib::Bond &b = path[i]->get_bond(*path[i + 1]);

				d += front_atom->crd().distance(b.atom2().crd());
				
				path[j]->set_max_linker_length(*path[i + 1], d);
				path[i + 1]->set_max_linker_length(*path[j], d);
				dbgmsg("last max_linker_length between " << *path[j] << " and " 
					<< *path[i + 1] << " = " << path[j]->get_max_linker_length(*path[i + 1]));
			}
			dbgmsg("TOTAL max_linker_length between " << *path[j] << " and " 
				<< *path[path.size() - 1] << " = " 
				<< path[j]->get_max_linker_length(*path[path.size() - 1]));
		}
	}

	
	void Segment::__set_branching_rules(const Segment::Paths &paths) {
		for (auto &kv : paths) {
			const Segment::Graph::Path &path = kv.second;
#ifndef NDEBUG
			dbgmsg("valid_path = ");
			for (auto &seg : path) dbgmsg(*seg);
#endif
			Segment &start = *path.back();
			Segment &goal = *path.front();
			Segment &start_next = **(path.end() - 2);
			Segment &goal_next = **(path.begin() + 1);
			const bool is_link_adjacent = __link_adjacent(path);
			for (auto it = path.begin(); it != path.end(); ++it) {
				Segment &current = **it;
				if (&current != &goal && goal.is_seed()) {
					Segment &next = **(it - 1);
					current.set_next(goal, next);
					goal.set_next(current, goal_next);
					if (is_link_adjacent)
						current.set_adjacent_seed_segments(goal);
					dbgmsg("current = " << current << " goal.is_seed() = " 
						<< boolalpha << goal.is_seed() << " next = " << next);
				}
				if (&current != &start && start.is_seed()) {
					Segment &prev = **(it + 1);
					current.set_next(start, prev);
					start.set_next(current, start_next);
					if (is_link_adjacent)
						current.set_adjacent_seed_segments(start);
					dbgmsg("current = " << current << " start.is_seed() = " 
						<< boolalpha << start.is_seed() << " prev = " << prev);
				}
			}
		}
	}
	

};
}
