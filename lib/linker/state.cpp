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
	
	//~ string State::pdb() const { 
		//~ stringstream ss;
		//~ for (auto &kv : __atom_crd) {
			//~ Molib::Atom a(*kv.first);
			//~ a.set_crd(kv.second);
			//~ ss << a;
		//~ } 
		//~ return ss.str();
	//~ }
	string State::pdb() const { 
		stringstream ss;
		for (int i = 0; i < get_segment().get_atoms().size(); ++i) {
			Molib::Atom a(get_segment().get_atom(i));
			a.set_crd(get_crd(i));
			ss << a;
		} 
		return ss.str();
	}
	
	bool State::clashes(const State &other, const Molib::Bond &excluded) const { // clashes between this and other state
		const double clash_coeff = 0.75;
		for (int i = 0; i < __crds.size(); ++i) {
			const Geom3D::Point &crd1 = get_crd(i);
			const Molib::Atom &a1 = __segment.get_atom(i);;
			if (&a1 == &excluded.atom1() || &a1 == &excluded.atom2()) { 
				dbgmsg("excluded state atom = " << a1.atom_number() 
					<< " is not checked for clashes");
				continue;
			}
			const double vdw1 = a1.radius();
			for (int j = 0; j < other.get_crds().size(); ++j) {
				const Geom3D::Point &crd2 = other.get_crd(j);
				const Molib::Atom &a2 = other.get_segment().get_atom(j);
				if (&a2 == &excluded.atom1() || &a2 == &excluded.atom2()) { 
					dbgmsg("excluded state atom = " << a2.atom_number() 
						<< " is not checked for clashes");
					continue;
				}
				const double vdw2 = a2.radius();
				if (crd1.distance_sq(crd2) < pow(clash_coeff * (vdw1 + vdw2), 2)) return true;
			}
		}
		return false;
	}
	
	//~ ostream& operator<< (ostream& stream, const State& s) {
		//~ stream << "State(address = " << &s <<", segment = " << s.__segment.get_seed_id() << ") " 
			//~ << " energy = " << setprecision(4) << fixed << s.__energy << " atom_crd =  " ;
		//~ for (auto &kv : s.__atom_crd) stream << kv.first->atom_number() << " -> " << kv.second << " ";
		//~ return stream;
	//~ }
	//~ 
	ostream& operator<< (ostream& stream, const State& s) {
		stream << "State(address = " << &s <<", segment = " << s.__segment.get_seed_id() << ") " 
			<< " energy = " << setprecision(4) << fixed << s.__energy << " atom_crd =  " ;
		for (int i = 0; i < s.get_crds().size(); ++i) {
			const Geom3D::Point &crd = s.get_crd(i);
			const Molib::Atom &a = s.__segment.get_atom(i);;
			stream << a.atom_number() << " -> " << crd << " ";
		}
		return stream;
	}
	
	ostream& operator<< (ostream& stream, const State::Vec& sv) {
		for (auto &state : sv) stream << "MEMBER STATE : " << state << endl;
		return stream;
	}
};
