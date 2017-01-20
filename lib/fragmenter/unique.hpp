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
			unique_ptr<BondGraph> graph;
			size_t seed_id;
		};
		typedef multimap<size_t, SeedData> USeeds;
		USeeds __unique_seeds;
		const string __seeds_file;
		void __read_seeds_file();
		bool __match(BondGraph&, USeeds::const_iterator, USeeds::const_iterator, size_t&) const;
		size_t __hash(const Atom::Set&) const;
		size_t __unique(const Atom::Set&);
	public:
		Unique(string seeds_file="") : __unique_seeds(USeeds()), __seeds_file(seeds_file) {
			if (__seeds_file != "") 
				__read_seeds_file();
		}
		~Unique() {
			if (__seeds_file != "") { // output to seeds_file if given
				stringstream ss;
				ss << __unique_seeds;
				inout::Inout::file_open_put_stream(__seeds_file, ss);
			}
		}
		size_t get_seed_id(const Atom::Set &a) { return __unique(a); }
		bool is_seed_unique( const Atom::Set &a) const;
		friend ostream& operator<<(ostream& os, const USeeds& useeds);
	};
}
#endif
