#ifndef DESIGN_DESIGN_H
#define DESIGN_DESIGN_H

#include "candock/molib/nrset.hpp"

namespace candock {

namespace design {

	class Design {
		molib::Molecules __designs;
		molib::Molecule  __original;
                molib::Unique &__existing;
	public:
		
		Design( const molib::Molecule &start, molib::Unique &existing );
		
		void functionalize_hydrogens_with_fragments(const molib::NRset& nr,
                                                            const double cutoff, const double clash_coeff,
                                                            const std::tuple<double, size_t, size_t, size_t>& lipinski_values
                                                           );
                static molib::Molecules functionalize_hydrogens_with_single_atoms( const molib::Molecule& original, const std::string& atom_type);
		void functionalize_extremes_with_single_atoms( const std::vector< std::string >& idatms );
		
		const molib::Molecules& get_internal_designs() const {
			return __designs;
		}
		
                void change_original_name( const std::string& name );
		
		const molib::Molecules& prepare_designs();
		const molib::Molecules& designs() const { return __designs; }
	};
}

}

#endif // DESIGN_DESIGN_H
