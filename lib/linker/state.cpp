#include "state.hpp"
#include "segment.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"

namespace Linker {
	
	State::Id State::idx = 0;

	State::Vec operator-(const State::Set& left, const State::Set& right) { 
		State::Vec ret; 
		for (auto &state : left) 
			if (!right.count(state)) 
				ret.push_back(state); 
		return ret; 
	}
	
	string State::pdb() const { 
		stringstream ss;
		for (auto &kv : __atom_crd) {
			Molib::Atom a(*kv.first);
			a.set_crd(kv.second);
			ss << a;
			//~ Atom &a = *const_cast<Atom*>(kv.first);
			//~ a.set_crd(kv.second);
			//~ ss << a;
		} 
		return ss.str();
	}
	
	bool State::clashes(const State &other, const Molib::Bond &excluded) const { // clashes between this and other state
		for (auto &kv1 : __atom_crd) {
			const Molib::Atom &a1 = *kv1.first;
			//~ if (&a1 == &excluded.first_atom() || &a1 == &excluded.second_atom()) { 
			if (&a1 == &excluded.atom1() || &a1 == &excluded.atom2()) { 
				dbgmsg("excluded state atom = " << a1.atom_number() 
					<< " is not checked for clashes");
				continue;
			}
			const Geom3D::Coordinate &c1 = kv1.second;
			const double vdw1 = a1.radius();
			for (auto &kv2 : other.get_atoms()) {
				const Molib::Atom &a2 = *kv2.first;
				//~ if (&a2 == &excluded.first_atom() || &a2 == &excluded.second_atom()) { 
				if (&a2 == &excluded.atom1() || &a2 == &excluded.atom2()) { 
					dbgmsg("excluded state atom = " << a2.atom_number() 
						<< " is not checked for clashes");
					continue;
				}
				const Geom3D::Coordinate &c2 = kv2.second;
				const double vdw2 = a2.radius();
				if (c1.distance_sq(c2) < pow(0.75 * (vdw1 + vdw2), 2)) return true;
			}
		}
		return false;
	}
	
	ostream& operator<< (ostream& stream, const State& s) {
		stream << "State(address = " << &s <<", segment = " << s.__segment.get_seed_id() << ") " 
			<< " energy = " << setprecision(4) << fixed << s.__energy << " atom_crd =  " ;
		for (auto &kv : s.__atom_crd) stream << kv.first->atom_number() << " -> " << kv.second << " ";
		return stream;
	}
	
	ostream& operator<< (ostream& stream, const State::Vec& sv) {
		for (auto &state : sv) stream << "MEMBER STATE : " << state << endl;
		return stream;
	}
};
