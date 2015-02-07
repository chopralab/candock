#include "segment.hpp"
//~ #include "linker.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"

namespace Molib {
	const Atom& Segment::adjacent_in_segment(const Atom &atom, 
		const Atom &forbidden) const { 
		//~ for (auto &bond : atom) {
			//~ auto &adj = bond.second_atom();
		for (auto &adj : atom) {
			if (&adj != &forbidden && has_atom(adj)) 
				return adj; 
		}
		throw Error("die : couldn't find adjacent in segment");
	}
};
