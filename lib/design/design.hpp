#ifndef DESIGN_DESIGN_H
#define DESIGN_DESIGN_H

#include "pdbreader/nrset.hpp"

namespace design {

	class Design {
		Molib::Molecules __designs;
		Molib::Molecule  __original;
	public:
		
		Design( const Molib::Molecule &start );
		
		void functionalize_hydrogens_with_fragments(const Molib::NRset& nr, const double cutoff, const double clash_coeff);
		void functionalize_hydrogens_with_single_atoms( const std::vector< std::string >& idatms );
		void functionalize_extremes_with_single_atoms( const std::vector< std::string >& idatms );
		
		const Molib::Molecules& get_internal_designs() const {
			return __designs;
		}
		
		const Molib::Molecules& prepare_designs( const std::string& seeds_file );
		const Molib::Molecules& designs() const { return __designs; }
	};
}

#endif // DESIGN_DESIGN_H