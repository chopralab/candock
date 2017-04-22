#ifndef UNIQUE_H
#include "helper/inout.hpp"
#include "fragmenter/fragmenter.hpp"
#include "molib/bond.hpp"
#include "molib/molecule.hpp"

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
			if (__seeds_file != "") { // output to seeds_file if given
				std::stringstream ss;
				ss << __unique_seeds;
				Inout::file_open_put_stream(__seeds_file, ss);
			}
		}
		size_t get_seed_id(const Atom::Set &a) { return __unique(a); }
		bool is_seed_unique( const Atom::Set &a) const;
		friend ostream& operator<<(ostream& os, const USeeds& useeds);
	};
}
#endif
