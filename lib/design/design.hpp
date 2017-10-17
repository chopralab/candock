#ifndef DESIGN_DESIGN_H
#define DESIGN_DESIGN_H

#include "molib/nrset.hpp"

namespace design {

	class Design {
		Molib::Molecules __designs;
		Molib::Molecule  __original;
                Molib::Unique &__existing;
	public:
		
		Design( const Molib::Molecule &start, Molib::Unique &existing );
		
		void functionalize_hydrogens_with_fragments(const Molib::NRset& nr,
                                                            const double cutoff, const double clash_coeff,
                                                            const std::tuple<double, size_t, size_t, size_t>& lipinski_values
                                                           );
		void functionalize_hydrogens_with_single_atoms( const std::vector< std::string >& idatms );
		void functionalize_extremes_with_single_atoms( const std::vector< std::string >& idatms );
		
		const Molib::Molecules& get_internal_designs() const {
			return __designs;
		}
		
                void change_original_name( const std::string& name );
		
		const Molib::Molecules& prepare_designs();
		const Molib::Molecules& designs() const { return __designs; }
	};
}

#endif // DESIGN_DESIGN_H
