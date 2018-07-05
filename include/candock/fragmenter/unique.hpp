#ifndef UNIQUE_H
#include "candock/helper/inout.hpp"
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/molib/bond.hpp"
#include "candock/molib/molecule.hpp"

namespace candock {

namespace Molib {
	class Atom;
	class Unique {
		struct SeedData {
			std::unique_ptr<BondGraph> graph;
			size_t seed_id;
		};
		typedef std::multimap<size_t, SeedData> USeeds;
		USeeds __unique_seeds;
		const std::string __seeds_file;
		void __read_seeds_file();
		bool __match(BondGraph&, USeeds::const_iterator, USeeds::const_iterator, size_t&) const;
		size_t __hash(const Atom::Set&) const;
		size_t __unique(const Atom::Set&);
	public:
		Unique(std::string seeds_file="") : __unique_seeds(USeeds()), __seeds_file(seeds_file) {
			if (__seeds_file != "") 
				__read_seeds_file();
		}

		~Unique() {
                        write_out();
		}

		size_t get_seed_id(const Atom::Set &a) { return __unique(a); }
		bool is_seed_unique( const Atom::Set &a) const;

                void write_out();

		friend ostream& operator<<(ostream& os, const USeeds& useeds);
	};
}

}

#endif
