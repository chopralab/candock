#ifndef UNIQUE_H
#define UNIQUE_H
#include "helper/help.hpp"
#include "helper/debug.hpp"
#include "helper/inout.hpp"
#include "graph/graph.hpp"
#include <tuple>
#include <functional>
#include "fragmenter/fragmenter.hpp"
#include "pdbreader/bond.hpp"
#include "pdbreader/molecule.hpp"


namespace Molib {
	class Atom;
	class Unique {
		struct SeedData {
			//~ unique_ptr<Glib::Graph<AtomTag>> graph;
			//~ unique_ptr<MolGraph> graph;
			//~ MolGraph graph;
			BondGraph graph;
			size_t seed_id;
		};
		typedef multimap<size_t, SeedData> USeeds;
		USeeds __unique_seeds;
		//~ size_t __seed_id;
		const string __seeds_file;
		void __read_seeds_file();
		//~ bool __match(Glib::Graph<AtomTag>&, USeeds::iterator, USeeds::iterator, size_t&);
		//~ bool __match(MolGraph&, USeeds::iterator, USeeds::iterator, size_t&);
		bool __match(BondGraph&, USeeds::iterator, USeeds::iterator, size_t&);
		size_t __hash(const AtomSet&);
		size_t __unique(const AtomSet&);
	public:
		//~ Unique(string seeds_file="") : __seeds_file(seeds_file), __seed_id(0), __unique_seeds(USeeds()) {
		Unique(string seeds_file="") : __seeds_file(seeds_file), 
			__unique_seeds(USeeds()) {
			if (__seeds_file != "") 
				__read_seeds_file();
		}
		//~ ~Unique() {
			//~ if (__seeds_file != "") { // output to seeds_file if given
				//~ stringstream ss;
				//~ for (auto &kv : __unique_seeds) {
					//~ ss << kv.second.seed_id << " " << kv.first << " " << kv.second.graph->get_smiles() << endl;
				//~ }
				//~ inout::Inout::file_open_put_stream(__seeds_file, ss);
			//~ }
		//~ }
		~Unique() {
			if (__seeds_file != "") { // output to seeds_file if given
				stringstream ss;
				ss << __unique_seeds;
				inout::Inout::file_open_put_stream(__seeds_file, ss);
			}
		}
		size_t get_seed_id(const AtomSet &a) { return __unique(a); }
		friend ostream& operator<<(ostream& os, const USeeds& useeds);
	};
	//~ AtomTags create_atom_tags(const help::smiles&);
}
#endif
