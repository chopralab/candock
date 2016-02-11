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
#include <map>
#include <algorithm>

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

	double Partial::compute_rmsd_ord(const Partial &other) const {
		Geom3D::Point::Vec crds1, crds2;
		map<const Molib::Atom*, pair<Geom3D::Point, Geom3D::Point>> crd_map;
		
		for (auto &pstate : this->get_states()) {
			for (int i = 0; i < pstate->get_segment().get_atoms().size(); ++i) {
				auto &a = pstate->get_segment().get_atom(i);
				crd_map[&a].first = pstate->get_crd(i);
			}
		}
		for (auto &pstate : other.get_states()) {
			for (int i = 0; i < pstate->get_segment().get_atoms().size(); ++i) {
				auto &a = pstate->get_segment().get_atom(i);
				crd_map[&a].second = pstate->get_crd(i);
			}
		}
		for (auto &kv : crd_map) {
			auto &crd1 = kv.second.first;
			auto &crd2 = kv.second.second;
			if (crd1 == Geom3D::Point() || crd2 == Geom3D::Point())
				throw Error("die : cannot compute rmsd between two partial conformations, which have different compositions");
			crds1.push_back(crd1);
			crds2.push_back(crd2);
		}
		return Geom3D::compute_rmsd(crds1, crds2);
	}

	Geom3D::Point Partial::compute_geometric_center() const { 
		return Geom3D::compute_geometric_center(this->get_ligand_crds()); 
	}

	void Partial::sort(Partial::Vec &v) {
		::sort(v.begin(), v.end(), Partial::comp());
	}

	Geom3D::Point::Vec Partial::get_ligand_crds() const { 
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
