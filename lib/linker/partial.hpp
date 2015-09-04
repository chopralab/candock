#ifndef CONFORMATION_H
#define CONFORMATION_H
#include "helper/debug.hpp"
#include "pdbreader/it.hpp"
#include "fragmenter/fragmenter.hpp"
#include "geom3d/coordinate.hpp"
#include "pdbreader/internal.hpp"
#include "pdbreader/grid.hpp"
#include <tuple>
#include <functional>
#include "segment.hpp"
#include "seed.hpp"
#include "helper/array2d.hpp"

namespace Linker {
	class State;

	class Partial {
	public:
		typedef vector<Partial> Vec;
		struct comp { bool operator()(const Partial &i, const Partial &j) 
			{ return i.get_energy() < j.get_energy(); } };

	private:
		State::Vec __states;
		double __energy;
		Geom3D::Point::Vec __crds;
	public:
		Partial() : __states(State::Vec()), __energy(0.0), __crds(Geom3D::Point::Vec()) {}
		Partial(const double energy) : __states(State::Vec()), __energy(energy), __crds(Geom3D::Point::Vec()) {}
		Partial(const State::Vec &states, const double energy, const Geom3D::Point::Vec &crds=Geom3D::Point::Vec())
			: __states(states), __energy(energy), __crds(crds) {}

		void add_state(State &state) { __states.push_back(&state); }

		State::Vec& get_states() { return __states; }
		const State::Vec& get_states() const { return __states; }

		void set_receptor_crds(const Geom3D::Point::Vec &crds) { __crds = crds; }
		Geom3D::Point::Vec& get_receptor_crds() { return __crds; }
		const Geom3D::Point::Vec& get_receptor_crds() const { return __crds; }

		void set_ligand_crds(const Geom3D::Point::Vec &crds); 
		Geom3D::Point::Vec get_ligand_crds() const;

		Molib::Atom::Vec get_ligand_atoms();

		void set_energy(const double energy) { __energy = energy; }
		double get_energy() const { return __energy; }

		int size() const { return __states.size(); }
		bool empty() const { return __states.empty(); }

		double compute_rmsd_ord(const Partial&) const;

		Geom3D::Point compute_geometric_center() const;
		
		static void sort(Partial::Vec &v);
		
		friend ostream& operator<<(ostream& os, const Partial &le);
		friend ostream& operator<<(ostream& os, const Vec &vec_le);
		

	};
}

#endif
