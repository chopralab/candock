#include "partial.hpp"
#include "state.hpp"
#include "score/score.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/bond.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"
#include "helper/array2d.hpp"
#include "graph/mcqd.hpp"
#include <queue>

using namespace std;

namespace Linker {

	ostream& operator<<(ostream& os, const Partial &le)	{
		os << "start link ++++++++++++++++++++++++++++++" << endl;
		os << "ENERGY = " << le.get_energy() << endl;
		for (auto &pstate : le.get_states())
			//~ os << *pstate << endl;
			os << pstate->pdb() << endl;
		os << "end link --------------------------------" << endl;
		return os;
	}

	ostream& operator<<(ostream& os, const Partial::Vec &vec_le)	{
		for (auto &le : vec_le) {
			os << le << endl;
		}			
		return os;
	}


	Geom3D::Point::Vec Partial::get_ligand_crds() { 
		Geom3D::Point::Vec crds;
		for (auto &pstate: __states) {
			for (auto &crd : pstate->get_crds()) {
				crds.push_back(crd);
			}
		}
		return crds; 
	}

	void Partial::set_ligand_crds(const Geom3D::Point::Vec &crds) { 
		int i = 0;
		for (auto &pstate: __states) {
			for (auto &crd : pstate->get_crds()) {
				crd = crds[i++];
			}
		}
	}

	Molib::Atom::Vec Partial::get_ligand_atoms() { 
		Molib::Atom::Vec atoms; 
		for (auto &pstate: __states) {
			for (auto &patom : pstate->get_segment().get_atoms()) {
				atoms.push_back(patom);
			}
		}
		return atoms;
	}

}
