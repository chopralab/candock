#ifndef ASSEMBLY_H
#define ASSEMBLY_H
#include "candock/geometry/geometry.hpp"
#include "candock/geometry/matrix.hpp"
#include "candock/molib/it.hpp"
#include "candock/molib/element.hpp"
#include "candock/molib/grid.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/residue.hpp"
#include "candock/molib/chain.hpp"
#include "candock/molib/model.hpp"

namespace candock {

namespace molib {
	class Chain;
	class Residue;
	class Atom;
	class Molecule;
	
	class Assembly : public template_map_container<Model, Assembly, Molecule> {
		int __number;
	std::string __name; // ASYMMETRIC UNIT OR BIOLOGICAL ASSEMBLY
	public:
		Assembly(int number, const std::string name="ASYMMETRIC UNIT") : __number(number), __name(name) {}
		Assembly(const Assembly &rhs) : __number(rhs.__number), __name(rhs.__name) { 
			for (auto &model : rhs) { 
				dbgmsg("Copy constructor : assembly");
				add(new Model(model));
			} 
		}
		
		void init_bio(const Assembly &asym, map<int, geometry::Matrix> &matrices, const set<char> &chains);
		
		void set_number(int number) { __number = number; }
		void set_name(const std::string &name) { __name = name; }
		int number() const { return __number; }
	std::string name() const { return __name; }
		Model& add(Model *m) { return this->aadd(m->number(), m, this); }

		void rotate(const geometry::Matrix &rota, const bool inverse=false);
		
		Atom::Vec get_atoms(const std::string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Assembly& erase_properties() { for (auto &model : *this) model.erase_properties(); return *this; }
		
		friend ostream& operator<< (ostream& stream, const Assembly& a);
	};
	

} // molib

}

#endif
	
