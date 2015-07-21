#ifndef STATE_H
#define STATE_H
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "geom3d/geom3d.hpp"
#include "pdbreader/bond.hpp"
#include "pdbreader/atom.hpp"
#include <tuple>
#include <functional>

namespace Molib {
	class Atom;
};

namespace Linker {
	class Segment;
	class State {
	public:
		typedef vector<State*> Vec;
		typedef vector<State> NVec;
		typedef set<State*> Set;
		typedef map<const State*, State*> Map;
		typedef vector<const State*> ConstVec;
		typedef pair<const State*, const State*> ConstPair;
		typedef int Id;
	private:
		static Id idx;
		const Segment &__segment;
		const Geom3D::Point::Vec __crds;
		double __energy;
		const Id __id;
		
	public:
		State(const Segment &segment, const Geom3D::Point::Vec crds, const double energy=0) : 
			__segment(segment), __crds(crds), __energy(energy), __id(idx++) { }
		void set_energy(const double energy) { __energy = energy; }
		double get_energy() const { return __energy; }			
		const Segment& get_segment() const { return __segment; }			
		const Geom3D::Point::Vec& get_crds() const { return __crds; }			
		const Geom3D::Point& get_crd(const int i) const { return __crds[i]; }
		bool clashes(const State &other, const Molib::Bond &excluded, const double clash_coeff) const; // clashes between this and other state
		string pdb() const;
		const Id get_id() const { return __id; }
		friend ostream& operator<< (ostream& stream, const State& s);
		friend ostream& operator<< (ostream& stream, const Vec& sv);
	};
	
	State::Vec operator-(const State::Set& left, const State::Set& right);
};
#endif
