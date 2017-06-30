#include "unique.hpp"
#include "molib/molecule.hpp"
#include "fragmenter/fragmenter.hpp"

namespace Molib {
	ostream& operator<<(ostream& os, const Unique::USeeds& useeds) {
		for (auto &kv : useeds) {
			os << kv.second.seed_id << " " << kv.first << " " 
                           << kv.second.graph->get_smiles() << endl;
		}
		return os;
	}

	void Unique::__read_seeds_file() {
		// seeds file is an index of unique seeds
		// each line is structured such as:
		// Seed_Id   hash(Chemical_Formula)  Mol_Graph
		// #seed_id   #hash(Car5 N1 ...)      Car#1_Car#2 #2_N2#3 ...
		vector<string> v_seeds;
		Inout::read_file(__seeds_file, v_seeds, Inout::no_panic);
		for (auto &s : v_seeds) {
			dbgmsg("reading seed file " << s);
			stringstream ss(s);
			size_t seed_id, hsh;
			ss >> seed_id >> hsh;
			dbgmsg(seed_id << " " << hsh);
			help::smiles edges;
			string vpair;
			while (ss >> vpair) {
				auto atom_props = help::ssplit(vpair, "_");
				edges.push_back(help::edge{atom_props[0], atom_props[1], ""});
			}
			__unique_seeds.insert(make_pair(hsh, 
				SeedData{unique_ptr<BondGraph>(new BondGraph(create_bonds(edges), true, false)), seed_id}));
		}
		dbgmsg("exiting read_seeds_file");
	}

	bool Unique::__match(BondGraph &g, USeeds::const_iterator it1, USeeds::const_iterator it2, size_t &si) const {
		dbgmsg("we are in __match");
		while (it1 != it2) {
			BondGraph &g2 = *(it1->second.graph);
			dbgmsg("before getting smiles");
			dbgmsg("g.size() == g2.size() " << boolalpha << (g.size() == g2.size()));
			dbgmsg("g.match(g2).size() " << g.match(g2)[0].first.size());
			if (g.isomorphic(g2)) {
				si = it1->second.seed_id;
				dbgmsg("finding equal seed number = " << si);
				return true;
			}
			it1++;
		}
		return false;
	}

	size_t Unique::__hash(const Atom::Set &atoms) const {
		map<string, int> chemical_formula;
		for (auto &a : atoms)
			chemical_formula[a->get_label()]++;
		stringstream ss;
		for (auto &kv : chemical_formula)
			ss << kv.first << " " << kv.second;
		std::hash<string> hash_fn;
		return hash_fn(ss.str());
	}

	size_t Unique::__unique(const Atom::Set &seed) {
		help::smiles edges;
		auto bonds = Molib::get_bonds_in(seed);
		for (auto &pbond : bonds) {
			const Molib::Bond &bond = *pbond;
			stringstream vertex1, vertex2;
			vertex1 << bond.atom1().get_label() << "#" << bond.atom1().atom_number();
			vertex2 << bond.atom2().get_label() << "#" << bond.atom2().atom_number();
			edges.push_back(help::edge{vertex1.str(), vertex2.str(), ""});
		}
		dbgmsg("before outputting edges");
		dbgmsg(edges);
		dbgmsg("before calculating hash");
		size_t hsh = __hash(seed);
		size_t si = 0;
		// Glib::Graph<AtomTag> g(create_atom_tags(s), true, false);
		// Atom::Graph g = create_graph(edges);
		dbgmsg("before creating bond graph");
		BondGraph g = create_graph(edges);
		dbgmsg(hsh);
		auto ret = __unique_seeds.equal_range(hsh);
#ifndef NDEBUG
		for (auto &kv : __unique_seeds) dbgmsg("hash = " << kv.first);
#endif
		dbgmsg("ret.first == ret.second " << boolalpha << (ret.first == ret.second));
		// seed's hash OR graph doesn't match any hash OR graph already in db, so add seed
		if (ret.first == ret.second || !__match(g, ret.first, ret.second, si)) { 
		//~ if (!__match(g, ret.first, ret.second, si)) { 
			// si = __seed_id++;
			si = __unique_seeds.size();
			dbgmsg("si = " << si);
			// __unique_seeds.insert(make_pair(hsh, 
				// SeedData{unique_ptr<Glib::Graph<AtomTag>>(new Glib::Graph<AtomTag>(create_atom_tags(s), true, false)), 
				// si}));
			__unique_seeds.insert(make_pair(hsh, 
				SeedData{unique_ptr<BondGraph>(new BondGraph(create_bonds(edges), true, false)), si}));
		}
		dbgmsg("si = " << si);
		return si;
	}
	bool Unique::is_seed_unique(const Atom::Set &seed) const {
		help::smiles edges;
		auto bonds = Molib::get_bonds_in(seed);
		for (auto &pbond : bonds) {
			const Molib::Bond &bond = *pbond;
			stringstream vertex1, vertex2;
			vertex1 << bond.atom1().get_label() << "#" << bond.atom1().atom_number();
			vertex2 << bond.atom2().get_label() << "#" << bond.atom2().atom_number();
			edges.push_back(help::edge{vertex1.str(), vertex2.str(), ""});
		}
		dbgmsg("before outputting edges");
		dbgmsg(edges);
		dbgmsg("before calculating hash");
		size_t hsh = __hash(seed);
		size_t si = 0;
		// Glib::Graph<AtomTag> g(create_atom_tags(s), true, false);
		// Atom::Graph g = create_graph(edges);
		dbgmsg("before creating bond graph");
		BondGraph g = create_graph(edges);
		dbgmsg(hsh);
		auto ret = __unique_seeds.equal_range(hsh);
#ifndef NDEBUG
		for (auto &kv : __unique_seeds) dbgmsg("hash = " << kv.first);
#endif
		dbgmsg("ret.first == ret.second " << boolalpha << (ret.first == ret.second));
		// seed's hash OR graph doesn't match any hash OR graph already in db, so add seed
		if (ret.first == ret.second || !__match(g, ret.first, ret.second, si)) { 
                        return true;
                }
		dbgmsg("si = " << si);
		return false;
	}

        void Unique::write_out() {
                if (__seeds_file != "") { // output to seeds_file if given
                        std::stringstream ss;
                        ss << __unique_seeds;
                        Inout::file_open_put_stream(__seeds_file, ss);
                }
        }

}
