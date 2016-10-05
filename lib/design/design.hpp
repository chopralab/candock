#ifndef DESIGN_DESIGN_H
#define DESIGN_DESIGN_H

#include "pdbreader/nrset.hpp"

namespace design {

	class Design {
		Molib::NRset __designs;
	public:
		static Molib::Molecules add_fragments_to_existing_molecule(const Molib::Molecule& start, const Molib::NRset& nr);
	};
}

#endif // DESIGN_DESIGN_H
