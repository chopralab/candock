/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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
                        write_out();
		}

		size_t get_seed_id(const Atom::Set &a) { return __unique(a); }
		bool is_seed_unique( const Atom::Set &a) const;

                void write_out();

		friend ostream& operator<<(ostream& os, const USeeds& useeds);
	};
}
#endif
