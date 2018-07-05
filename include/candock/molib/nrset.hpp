#ifndef NRSET_H
#define NRSET_H
#include "candock/molib/it.hpp"
#include "candock/molib/residue.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/molib/molecules.hpp"

#include <random>

namespace candock {

namespace Molib {
	
	class NRset : public template_map_container<Molecules, NRset, NRset> {
	public:
		Molecules& add(Molecules *m) { return this->aadd(m, this); }
		NRset& erase_properties() { for (auto &molecules : *this) molecules.erase_properties(); return *this; }
		Molecule::Vec get_molecules(const Residue::res_type &rest) const;
		Atom::Vec get_atoms(const std::string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		void jiggle(std::mt19937 rng);
		friend ostream& operator<< (ostream& stream, const NRset& m);
	};

} // Molib

}

#endif
