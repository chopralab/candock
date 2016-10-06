#ifndef DESIGN_DESIGN_H
#define DESIGN_DESIGN_H

#include "pdbreader/nrset.hpp"

namespace design {

	class Design {
		Molib::Molecules __designs;
		Molib::Molecule  __original;
	public:
		
		Design( const Molib::Molecule &start );
		
		void add_fragments_to_existing_molecule(const Molib::NRset& nr);
		
		const Molib::Molecules& get_internal_designs() const {
			return __designs;
		}
		
		const Molib::Molecules& get_prepared_designs();
	};
}

#endif // DESIGN_DESIGN_H
