#include "unique.hpp"
#include "pdbreader/molecule.hpp"
#include "fragmenter/fragmenter.hpp"

namespace Molib {
	ostream& operator<<(ostream& os, const Unique::USeeds& useeds) {
		for (auto &kv : useeds) {
			os << kv.second.seed_id << " " << kv.first << " " 
			<< kv.second.graph.get_smiles() << endl;
		}
	}
	//~ void Unique::__read_seeds_file() {
		//~ // seeds file is an index of unique seeds
		//~ // each line is structured such as:
		//~ // Seed_Id   hash(Chemical_Formula)  Mol_Graph
		//~ // #seed_id   #hash(Car5 N1 ...)      Car#1_Car#2 #2_N2#3 ...
		//~ vector<string> v_seeds;
		//~ inout::Inout::read_file(__seeds_file, v_seeds, inout::Inout::no_panic);
		//~ for (auto &s : v_seeds) {
			//~ dbgmsg("reading seed file " << s);
			//~ stringstream ss(s);
			//~ size_t seed_id, hsh;
			//~ ss >> seed_id >> hsh;
			//~ dbgmsg(seed_id << " " << hsh);
			//help::smiles bonds;
			//~ help::smiles edges;
			//~ string vpair;
			//~ while (ss >> vpair) {
				//bonds.push_back(pair<string, string>(help::ssplit(vpair, "_")[0], help::ssplit(vpair, "_")[1]));
				//edges.push_back(pair<string, string>(help::ssplit(vpair, "_")[0], help::ssplit(vpair, "_")[1]));
				//~ auto atom_props = help::ssplit(vpair, "_");
				//~ edges.push_back(help::edge{atom_props[0], atom_props[1], ""});
				//~ dbgmsg(vpair);
			//~ }
			//__unique_seeds.insert(make_pair(hsh, 
			//	SeedData{unique_ptr<Glib::Graph<AtomTag>>(new Glib::Graph<AtomTag>(create_atom_tags(bonds), true, false)), 
			//__unique_seeds.insert(make_pair(hsh, SeedData{create_graph(edges), seed_id}));
			//~ __unique_seeds.insert({hsh, SeedData{create_graph(edges), seed_id}});
		//~ }
		//~ dbgmsg("exiting read_seeds_file");
	//~ }
	void Unique::__read_seeds_file() {
		// seeds file is an index of unique seeds
		// each line is structured such as:
		// Seed_Id   hash(Chemical_Formula)  Mol_Graph
		// #seed_id   #hash(Car5 N1 ...)      Car#1_Car#2 #2_N2#3 ...
		vector<string> v_seeds;
		inout::Inout::read_file(__seeds_file, v_seeds, inout::Inout::no_panic);
		for (auto &s : v_seeds) {
			dbgmsg("reading seed file " << s);
			stringstream ss(s);
			size_t seed_id, hsh;
			ss >> seed_id >> hsh;
			dbgmsg(seed_id << " " << hsh);
			//~ help::smiles bonds;
			help::smiles edges;
			string vpair;
			while (ss >> vpair) {
				//~ bonds.push_back(pair<string, string>(help::ssplit(vpair, "_")[0], help::ssplit(vpair, "_")[1]));
				//~ edges.push_back(pair<string, string>(help::ssplit(vpair, "_")[0], help::ssplit(vpair, "_")[1]));
				auto atom_props = help::ssplit(vpair, "_");
				edges.push_back(help::edge{atom_props[0], atom_props[1], ""});
				//~ dbgmsg(vpair);
			}
			//~ __unique_seeds.insert(make_pair(hsh, 
				//~ SeedData{unique_ptr<Glib::Graph<AtomTag>>(new Glib::Graph<AtomTag>(create_atom_tags(bonds), true, false)), 
			//~ __unique_seeds.insert(make_pair(hsh, SeedData{create_graph(edges), seed_id}));
			//~ __unique_seeds.insert({hsh, SeedData{create_graph(edges), seed_id}});
			dbgmsg("before create_mol_graph edges = " << endl << edges);
			__unique_seeds.insert({hsh, SeedData{create_mol_graph(edges), seed_id}});
			//~ create_mol_graph(edges);
			//~ __unique_seeds.insert({hsh, 
				//~ SeedData{unique_ptr<MolGraph>(new MolGraph(create_mol_graph(edges))), seed_id}});
		}
		throw Error("exit after read seeds file");
		dbgmsg("exiting read_seeds_file");
	}
	//~ bool Unique::__match(Glib::Graph<AtomTag> &g, USeeds::iterator it1, USeeds::iterator it2, size_t &si) {
	//~ bool Unique::__match(MolGraph &g, USeeds::iterator it1, USeeds::iterator it2, size_t &si) {
	//~ bool Unique::__match(BondGraph &g, USeeds::iterator it1, USeeds::iterator it2, size_t &si) {
		//~ dbgmsg("we are in __match");
		//~ while (it1 != it2) {
			//~ BondGraph &g2 = it1->second.graph;
			//~ dbgmsg("before getting smiles");
			//~ dbgmsg(g.get_smiles());
			//~ dbgmsg(g2.get_smiles());
			//~ dbgmsg("g.size() == g2.size() " << boolalpha << (g.size() == g2.size()));
			//~ dbgmsg("g.match(g2).size() " << g.match(g2)[0].first.size());
			//~ if (g.isomorphic(g2)) {
				//~ si = it1->second.seed_id;
				//~ dbgmsg("finding equal seed number = " << si);
				//~ return true;
			//~ }
			//~ it1++;
		//~ }
		//~ return false;
	//~ }
	bool Unique::__match(MolGraph &g, USeeds::iterator it1, USeeds::iterator it2, size_t &si) {
		dbgmsg("we are in __match");
		while (it1 != it2) {
			//~ Glib::Graph<AtomTag> &g2 = *it1->second.graph;
			MolGraph &g2 = it1->second.graph;
			//~ BondGraph &g2 = it1->second.graph;
			dbgmsg("before getting smiles");
			dbgmsg(g.get_smiles());
			dbgmsg(g2.get_smiles());
			dbgmsg("g.size() == g2.size() " << boolalpha << (g.size() == g2.size()));
			dbgmsg("g.match(g2).size() " << g.match(g2)[0].first.size());
			if (g.isomorphic(g2)) {
				si = it1->second.seed_id;
				//~ size_t si = it1->second.seed_id;
				dbgmsg("finding equal seed number = " << si);
				return true;
			}
			it1++;
		}
		return false;
	}
	size_t Unique::__hash(const AtomSet &atoms) {
		map<string, int> chemical_formula;
		for (auto &a : atoms)
			chemical_formula[a->get_label()]++;
		stringstream ss;
		for (auto &kv : chemical_formula)
			ss << kv.first << " " << kv.second;
		std::hash<string> hash_fn;
		return hash_fn(ss.str());
	}
	//~ size_t Unique::__unique(const AtomSet &seed) {
		//~ help::smiles edges;
		//~ for (auto &atom : seed) {
			//~ const Molib::Atom &a = *atom;
			//~ for (auto &adj_a : a) {
				//~ if (a.atom_number() < adj_a.atom_number() && seed.find(&adj_a) != seed.end()) {
					//~ stringstream vertex1, vertex2;
					//~ vertex1 << a.get_label() << "#" << a.atom_number();
					//~ vertex2 << adj_a.get_label() << "#" << adj_a.atom_number();
					//~ edges.push_back(help::edge{vertex1.str(), vertex2.str(), ""});
				//~ }
			//~ }
		//~ }
		//~ dbgmsg("before outputting edges");
		//~ dbgmsg(edges);
		//~ dbgmsg("before calculating hash");
		//~ size_t hsh = __hash(seed);
		//~ size_t si = 0;
		//~ // Glib::Graph<AtomTag> g(create_atom_tags(s), true, false);
		//~ // MolGraph g = create_graph(edges);
		//~ dbgmsg("before creating bond graph");
		//~ BondGraph g = create_graph(edges);
		//~ dbgmsg(hsh);
		//~ auto ret = __unique_seeds.equal_range(hsh);
//~ #ifndef NDEBUG
		//~ for (auto &kv : __unique_seeds) dbgmsg("hash = " << kv.first);
//~ #endif
		//~ dbgmsg("ret.first == ret.second " << boolalpha << (ret.first == ret.second));
		//~ // seed's hash OR graph doesn't match any hash OR graph already in db, so add seed
		//~ if (ret.first == ret.second || !__match(g, ret.first, ret.second, si)) { 
			//~ // si = __seed_id++;
			//~ si = __unique_seeds.size();
			//~ dbgmsg(si);
			//~ // __unique_seeds.insert(make_pair(hsh, 
				//~ // SeedData{unique_ptr<Glib::Graph<AtomTag>>(new Glib::Graph<AtomTag>(create_atom_tags(s), true, false)), 
				//~ // si}));
			//~ __unique_seeds.insert(make_pair(hsh, 
				//~ SeedData{create_graph(edges), si}));
		//~ }
		//~ dbgmsg(si);
		//~ return si;
	//~ }
	size_t Unique::__unique(const AtomSet &seed) {
		//~ help::smiles edges;
		//~ for (auto &patom : seed) {
			//~ const Molib::Atom &atom = *patom;
			//~ for (auto &adj_a : a) {
				//~ if (a.atom_number() < adj_a.atom_number() && seed.find(&adj_a) != seed.end()) {
					//~ stringstream vertex1, vertex2;
					//~ vertex1 << a.get_label() << "#" << a.atom_number();
					//~ vertex2 << adj_a.get_label() << "#" << adj_a.atom_number();
					//~ edges.push_back(help::edge{vertex1.str(), vertex2.str(), ""});
				//~ }
			//~ }
		//~ }
		//~ dbgmsg("before outputting edges");
		//~ dbgmsg(edges);
//~ #ifndef NDEBUG
		//~ for (auto &kv : s)
			//~ dbgmsg("vertex1 = " << kv.first << " vertex2 = " << kv.second);
//~ #endif
		dbgmsg("before calculating hash");
		size_t hsh = __hash(seed);
		size_t si = 0;
		//~ Glib::Graph<AtomTag> g(create_atom_tags(s), true, false);
		//~ MolGraph g = create_graph(edges);
		dbgmsg("before creating atom graph");
		//~ BondGraph g = create_graph(edges);
		MolGraph g = create_graph(seed);
		dbgmsg(hsh);
		auto ret = __unique_seeds.equal_range(hsh);
#ifndef NDEBUG
		for (auto &kv : __unique_seeds) dbgmsg("hash = " << kv.first);
#endif
		dbgmsg("ret.first == ret.second " << boolalpha << (ret.first == ret.second));
		// seed's hash OR graph doesn't match any hash OR graph already in db, so add seed
		if (ret.first == ret.second || !__match(g, ret.first, ret.second, si)) { 
			//~ si = __seed_id++;
			si = __unique_seeds.size();
			dbgmsg(si);
			//~ __unique_seeds.insert(make_pair(hsh, 
				//~ SeedData{unique_ptr<Glib::Graph<AtomTag>>(new Glib::Graph<AtomTag>(create_atom_tags(s), true, false)), 
				//~ si}));
			//~ __unique_seeds.insert(make_pair(hsh, 
				//~ SeedData{create_graph(edges), si}));
			__unique_seeds.insert({hsh,	SeedData{g, si}});
		}
		dbgmsg(si);
		return si;
	}
}
