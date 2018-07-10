#ifndef CONFORMATION_H
#define CONFORMATION_H
#include "candock/helper/debug.hpp"
#include "candock/molib/it.hpp"
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/geometry/coordinate.hpp"
#include "candock/molib/internal.hpp"
#include "candock/molib/grid.hpp"
#include "candock/linker/segment.hpp"
#include "candock/linker/seed.hpp"
#include "candock/helper/array2d.hpp"

#include <tuple>
#include <functional>

namespace candock {

namespace linker {
	class State;

	class Partial {
	public:
		typedef std::vector<Partial> Vec;
		struct comp { bool operator()(const Partial &i, const Partial &j) 
			{ return i.get_energy() < j.get_energy(); } };

	private:
		State::Vec __states;
		double __energy;
		geometry::Point::Vec __crds;
	public:
		Partial() : __states(State::Vec()), __energy(0.0), __crds(geometry::Point::Vec()) {}
		Partial(const double energy) : __states(State::Vec()), __energy(energy), __crds(geometry::Point::Vec()) {}
		Partial(const State::Vec &states, const double energy, const geometry::Point::Vec &crds=geometry::Point::Vec())
			: __states(states), __energy(energy), __crds(crds) {}

		void add_state(State &state) { __states.push_back(&state); }

		State::Vec& get_states() { return __states; }
		const State::Vec& get_states() const { return __states; }

		void set_receptor_crds(const geometry::Point::Vec &crds) { __crds = crds; }
		geometry::Point::Vec& get_receptor_crds() { return __crds; }
		const geometry::Point::Vec& get_receptor_crds() const { return __crds; }

		void set_ligand_crds(const geometry::Point::Vec &crds); 
		geometry::Point::Vec get_ligand_crds() const;

		molib::Atom::Vec get_ligand_atoms();

		void set_energy(const double energy) { __energy = energy; }
		double get_energy() const { return __energy; }

		int size() const { return __states.size(); }
		bool empty() const { return __states.empty(); }

		double compute_rmsd_ord(const Partial&) const;

		geometry::Point compute_geometric_center() const;
		
		static void sort(Partial::Vec &v);
		
		friend std::ostream& operator<<(std::ostream& os, const Partial &le);
		friend std::ostream& operator<<(std::ostream& os, const Vec &vec_le);
		

	};
}

}

#endif
