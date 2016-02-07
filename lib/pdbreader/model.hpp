#ifndef MODEL_H
#define MODEL_H
#include "geom3d/geom3d.hpp"
#include "geom3d/matrix.hpp"
#include "fragmenter/fragmenter.hpp"
#include "it.hpp"
#include "element.hpp"
#include "grid.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"

using namespace std;

namespace Molib {
	class Chain;
	class Residue;
	class Atom;
	class Assembly;
	
	class Model : public template_map_container<Chain, Model, Assembly, char> {
		int __number;
		map<pair<int, Residue::res_tuple2>, map<Residue::res_tuple2, Residue::res_tuple2>> __remarks;
		Fragmenter::Fragment::Vec __rigid;
	public:
		Model(int number) : __number(number) {}
		Model(const Model &rhs) : __number(rhs.__number), __remarks(rhs.__remarks) { 
			for (auto &chain : rhs) { 
				dbgmsg("Copy constructor : model");
				add(new Chain(chain)); 
			} 
		}
		Fragmenter::Fragment::Vec& get_rigid() { return __rigid; }
		const Fragmenter::Fragment::Vec& get_rigid() const { return __rigid; }
		void set_rigid(const Fragmenter::Fragment::Vec &rigid) { __rigid = rigid; }
		void init_bio(const Model &model_asym, const Geom3D::Matrix &matrix, const set<char> &chains);
		void rotate(const Geom3D::Matrix &rota, const bool inverse=false);
		Chain& add(Chain *chain) { return this->aadd(chain->chain_id(), chain, this); }
		void set_remark(const int remark_number, const Residue::res_tuple2 &ligand, pair<const Residue::res_tuple2&, const Residue::res_tuple2&> rpair) { 
				this->__remarks[make_pair(remark_number, ligand)]
					.insert(make_pair(rpair.first, rpair.second)); 
		}
		Model& regenerate_bonds(const Model&);
		void set_number(const int &i) { __number = i; }
		Chain& chain(const char chain_id) const { return this->element(chain_id); }
		bool has_chain(const char chain_id) const { return this->has_element(chain_id); }
		int number() const { return __number; }
		
		bool remarks(const int remark_number, const Residue::res_tuple2 ligand, map<Residue::res_tuple2, Residue::res_tuple2> &r) { 
			auto i = __remarks.find(make_pair(remark_number, ligand));
			if (i == __remarks.end()) return false;
			r = i->second;
			return true;
		}
		
		Atom::Vec get_atoms(const string &chain_ids="", const Residue::res_type &rest=Residue::res_type::notassigned) const;
		Model& erase_properties() { for (auto &chain : *this) chain.erase_properties(); return *this; }
		friend ostream& operator<< (ostream& stream, const Model& m);
	};

} // Molib
#endif
	
