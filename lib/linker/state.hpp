#ifndef STATE_H
#define STATE_H
#include "helper/debug.hpp"
#include "helper/help.hpp"
#include "geom3d/coordinate.hpp"
#include "pdbreader/bond.hpp"
#include <tuple>
#include <functional>

namespace Molib {
	class Atom;
	class Segment;
	typedef map<const Atom*, Geom3D::Coordinate> AtomToCrd;
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
		static Id idx = 0;
		const Segment &__segment;
		const AtomToCrd __atom_crd;
		double __energy;
		Id __id;
		
	public:
		State(const Segment &segment, const AtomToCrd atom_crd, const double energy=0) : 
			__segment(segment), __atom_crd(atom_crd), __energy(energy), __id(idx++) { }
		void set_energy(const double energy) { __energy = energy; }
		double get_energy() const { return __energy; }			
		const Segment& get_segment() const { return __segment; }			
		const AtomToCrd& get_atoms() const { return __atom_crd; }			
		Geom3D::Coordinate get_atom_crd(const Atom &atom) const { return __atom_crd.at(&atom); }
		bool has_atom(const Atom &atom) const { return __atom_crd.count(&atom); }
		bool clashes(const State &other, const Bond &excluded) const; // clashes between this and other state
		string pdb() const;
		friend ostream& operator<< (ostream& stream, const State& s);
		friend ostream& operator<< (ostream& stream, const Vec& sv);
	};
	
	State::Vec operator-(const State::Set& left, const State::Set& right) const;
};
#endif
