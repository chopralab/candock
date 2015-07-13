#include "helper/inout.hpp"
#include "helper/debug.hpp"
#include "helper/benchmark.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "geom3d/geom3d.hpp"
#include "graph/mcqd.hpp"
#include "helper/array2d.hpp"
#include "gpoints.hpp"
#include "conformations.hpp"
#include <iostream>
#include <exception>
using namespace std;


namespace Docker {

	ostream& operator<<(ostream& os, const Conformations::Conf &conf) {
		os << "MODEL" << endl;
		for (int i = 0; i < conf.get_points().size(); ++i) {
			auto &gpoint0 = conf.get_point(i);
			auto &atom = conf.get_atom(i);
			atom.set_crd(gpoint0.crd());
			os << atom;
		}
		os << "ENDMDL" << endl;
		return os;
	}

	ostream& operator<<(ostream& os, const Conformations &conformations) {
		for (auto &conf : conformations.__conf_vec) {
			os << conf << endl;
		}
		return os;
	}

	double Conformations::Conf::compute_rmsd(const Conformations::Conf &other) {
		
		double sum_sq(0);
		for (int i = 0; i < this->__points.size(); ++i) {
			Geom3D::Point &crd1 = this->__points[i]->crd();
			Geom3D::Point &crd2 = other.__points[i]->crd();
			sum_sq += crd1.distance_sq(crd2);
		}
		return sqrt(sum_sq / this->__points.size());
	}

	Conformations::Conformations(const Molib::Molecule &seed, Gpoints &gpoints, 
			const double &grid_spacing) {
		try {
			__init_conformations(seed, gpoints, grid_spacing);
		} catch(...) {
			dbgmsg("FAILURE: constructor of Conformations failed ... cleaning resources...");
			throw;
		}
	}
	
	void Conformations::__init_conformations(const Molib::Molecule &seed, Gpoints &gpoints, 
		const double &grid_spacing) {
	
		Benchmark::reset();

		// order seed atoms according to their atom numbers
		__ordered_atoms = seed.get_atoms();
		sort(__ordered_atoms.begin(), __ordered_atoms.end(), 
			[](const Molib::Atom *i, const Molib::Atom *j) { return i->atom_number() < j->atom_number(); }); 

		// create a mapping between atom pointers and their positions in vec
		map<Molib::Atom*, int> atom_to_i;
		for (int i = 0; i < __ordered_atoms.size(); ++i) atom_to_i[__ordered_atoms[i]] = i;

		// center the grid points around the center point
		Docker::Gpoints::Gpoint cp = gpoints.get_center_point();
		dbgmsg("center point = " << cp.ijk());
		for (auto &point : gpoints.get_gridpoints0()) {
			dbgmsg("point = " << point.ijk());
			point.ijk() = point.ijk() - cp.ijk();
			dbgmsg("centered point = " << point.ijk());
		}
		
		// create grid
		Docker::Gpoints::PGpointVec pgvec;
		for (auto &point : gpoints.get_gridpoints0()) {
			pgvec.push_back(&point);
		 }
		Grid<Docker::Gpoints::Gpoint> grid(pgvec);
		dbgmsg("after creating grid");

		// get maximum distance from atom A in seed to any atom B
		Molib::Atom &atomA = seed.get_center_atom();
		Molib::AtomVec atoms;
		dbgmsg("atomA = " << atomA);
		
		// get seed atoms without atomA
		for (auto &patom : __ordered_atoms) {
			if (patom != &atomA) atoms.push_back(patom);
		}
		
		// get point that is at the center of the gpoints (0,0,0)
		Docker::Gpoints::Gpoint &pointA = gpoints.get_center_point();

		// generate product graph vertices
		vector<pair<Molib::Atom*, Docker::Gpoints::Gpoint*>> vertices;
		vertices.push_back(make_pair(&atomA, &pointA));
		const double tol = grid_spacing / 3;

		for (auto &patomB : atoms) {
			for (auto &pointB : grid.get_neighbors_within_tolerance(pointA, atomA.crd().distance(patomB->crd()), tol)) {
				vertices.push_back(make_pair(patomB, pointB));
			}
		}
		dbgmsg("number of vertices = " << vertices.size());

		// set the adjacency matrix (generate product graph edges)
		Array2d<bool> conn(vertices.size(), vertices.size());
		for (int i = 0; i < vertices.size(); ++i) {
			auto &v1 = vertices[i];
			for (int j = i + 1; j < vertices.size(); ++j) {
				auto &v2 = vertices[j];
				if (v1.first != v2.first && v1.second != v2.second 
					&& fabs(v1.first->crd().distance(v2.first->crd()) - v1.second->crd().distance(v2.second->crd())) < tol) {
					conn.data[i][j] = conn.data[j][i] = true;
				}
			}
		}
		
		// find all max cliques of size equal num. seed atoms
		Maxclique m(conn.data, vertices.size());
		vector<vector<int>> qmaxes = m.mcq(__ordered_atoms.size());

		// save conformations
		for (int c = 0; c < qmaxes.size(); ++c) {
			auto &qmax = qmaxes[c];
			// if max clique contains all seed atoms...
			dbgmsg("qmax[" << c << "].size = " << qmax.size());
			// go over vertices of each max clique
			Docker::Gpoints::PGpointVec points(qmax.size()); // init size of vec!
			for (auto &i : qmax) {
				Molib::Atom *patom = vertices[i].first;
				Docker::Gpoints::Gpoint *gpoint = vertices[i].second;
				Docker::Gpoints::IJK ijk = gpoint->ijk();
				__conf_map[ijk.i][ijk.j][ijk.k].push_back(c);
				points[atom_to_i[patom]] = gpoint;
			}
			//~ __conf_vec.push_back(Conformations::Conf(__atom_matches, __ordered_atoms, points));
			__conf_vec.push_back(Conformations::Conf(__ordered_atoms, points));
		}
		
		// pre-compute rmsd between ALL conformations
		__rmsd.init(__conf_vec.size(), __conf_vec.size());
		for (int i = 0; i < __conf_vec.size(); ++i) {
			for (int j = i + 1; j < __conf_vec.size(); ++j) {
				__rmsd.data[i][j] = __rmsd.data[j][i] = __conf_vec[i].compute_rmsd(__conf_vec[j]);
			}
		}
		cout << "time to find " << __conf_vec.size() << " conformations of seed " 
			<< seed.name() << " took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
	}
	
};