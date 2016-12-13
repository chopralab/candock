#ifndef ASSEMBLY_H
#define ASSEMBLY_H
#include "geom3d/geom3d.hpp"
#include "geom3d/matrix.hpp"
#include "it.hpp"
#include "element.hpp"
#include "grid.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
#include "model.hpp"

using namespace std;

namespace Molib {
	class Chain;
	class Residue;
	class Atom;
	class Molecule;
	
	class Assembly : public template_map_container<Model, Assembly, Molecule> {
		int __number;
		string __name; // ASYMMETRIC UNIT OR BIOLOGICAL ASSEMBLY
	public:
		Assembly(int number, const string name="ASYMMETRIC UNIT") : __number(number), __name(name) {}
		Assembly(const Assembly &rhs) : __number(rhs.__number), __name(rhs.__name) { 
			for (auto &model : rhs) { 
				dbgmsg("Copy constructor : assembly");
				add(new Model(model));
			} 
		}
		
		void init_bio(const Assembly &asym, map<int, Geom3D::Matrix> &matrices, const set<char> &chains);
		
		void set_number(int number) { __number = number; }
		void set_name(const string &name) { __name = name; }
		int number() const { return __number; }
		string name() const { return __name; }
		Model& add(Model *m) { return this->aadd(m->number(), m, this); }

		void rotate(const Geom3D::Matrix &rota, const bool inverse=false);
		
		Atom::Vec get_atoms(const string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned, const int model_number=-1) const;
		Assembly& erase_properties() { for (auto &model : *this) model.erase_properties(); return *this; }
		
		friend ostream& operator<< (ostream& stream, const Assembly& a);
	};
	

} // Molib
#endif
	
