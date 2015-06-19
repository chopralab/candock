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

ostream& operator<<(ostream& os, const Docker::OneConformation &conf) {
	os << "MODEL" << endl;
	for (auto &pair : conf) {
		Molib::Atom &atom = *pair.first;
		Docker::Gpoint &gpoint0 = *pair.second;
		atom.set_crd(gpoint0.crd());
		os << atom;
	}
	os << "ENDMDL" << endl;
	return os;
}


namespace Docker {

	ostream& operator<<(ostream& os, const Conformations &conformations) {
		for (auto &conf : conformations.__conf) {
			os << conf << endl;
		}
		return os;
	}

	Conformations::Conformations(const Molib::Molecule &molecule, Gpoints &gpoints, 
			const double &grid_spacing) {
		try {
			__init_conformations(molecule, gpoints, grid_spacing);
		} catch(...) {
			dbgmsg("FAILURE: constructor of Conformations failed ... cleaning resources...");
			throw;
		}
	}
	
	void Conformations::__init_conformations(const Molib::Molecule &seed, Gpoints &gpoints, 
		const double &grid_spacing) {
	
		Benchmark::reset();
			
		Docker::Gpoint cp = gpoints.get_center_point();
		dbgmsg("center point = " << cp.ijk());
		for (auto &point : gpoints.get_gridpoints0()) {
			dbgmsg("point = " << point.ijk());
			point.ijk() = point.ijk() - cp.ijk();
			dbgmsg("centered point = " << point.ijk());
		}
		// create grid
		vector<Docker::Gpoint*> pgpvec;
		for (auto &point : gpoints.get_gridpoints0()) {
			pgpvec.push_back(&point);
		 }
		Grid<Docker::Gpoint> grid(pgpvec);
		dbgmsg("after creating grid");
		// get maximum distance from atom A in seed to any atom B
		const double tol = grid_spacing / 3;
		Molib::Atom &atomA = seed.get_center_atom();
		Molib::AtomVec atoms;
		dbgmsg("atomA = " << atomA);
		// get seed atoms without atomA
		for (auto &patom : seed.get_atoms()) {
			if (patom != &atomA) atoms.push_back(patom);
		}

		// get point that is at the center of the gpoints (0,0,0)
		Docker::Gpoint &pointA = gpoints.get_center_point();

		vector<pair<Molib::Atom*, Docker::Gpoint*>> vertices;
		vertices.push_back(make_pair(&atomA, &pointA));
		// generate product graph vertices
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
		vector<vector<int>> qmaxes = m.mcq(seed.get_atoms().size());

		// save conformations
		for (int c = 0; c < qmaxes.size(); ++c) {
			auto &qmax = qmaxes[c];
			// if max clique contains all seed atoms...
			dbgmsg("qmax[" << c << "].size = " << qmax.size());
			OneConformation one_conf;
			// go over vertices of each max clique
			for (auto &i : qmax) {
				Molib::Atom *patom = vertices[i].first;
				Docker::Gpoint *gpoint = vertices[i].second;
				Docker::IJK ijk = gpoint->ijk();
				__confmap[ijk.i][ijk.j][ijk.k].push_back(c);
				one_conf.push_back(make_pair(patom, gpoint));
			}
			__conf.push_back(one_conf);
		}
		cout << "time to find " << __conf.size() << " conformations of seed " 
			<< seed.name() << " took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
	}
	
};
