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
		
		void functionalize_hydrogens_with_fragments(const Molib::NRset& nr, const double cutoff, const double clash_coeff);
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
