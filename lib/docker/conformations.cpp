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

	Conformations::Conformations(const Molib::Molecule &seed, const Gpoints &gpoints, 
		const double &grid_spacing, const int min_num_conf) {
		try {
			__init_conformations(seed, gpoints, grid_spacing, min_num_conf);
		} catch(...) {
			dbgmsg("FAILURE: constructor of Conformations failed ... cleaning resources...");
			throw;
		}
	}
	
	void Conformations::__init_conformations(const Molib::Molecule &seed, const Gpoints &gpoints, 
		const double &grid_spacing, const int min_num_conf) {
	
		Benchmark::reset();

		// order seed atoms according to their atom numbers
		__ordered_atoms = seed.get_atoms();
		sort(__ordered_atoms.begin(), __ordered_atoms.end(), 
			[](const Molib::Atom *i, const Molib::Atom *j) { return i->atom_number() < j->atom_number(); }); 

		// create a mapping between atom pointers and their positions in vec
		map<Molib::Atom*, int> atom_to_i;
		for (int i = 0; i < __ordered_atoms.size(); ++i) atom_to_i[__ordered_atoms[i]] = i;

		// create grid
		Gpoints::PGpointVec pgvec;
		for (auto &point : gpoints.get_gridpoints0()) {
			pgvec.push_back(const_cast<Gpoints::Gpoint*>(&point));
		 }
		Grid<Gpoints::Gpoint> grid(pgvec);
		dbgmsg("after creating grid");

		// get maximum distance from atom A in seed to any atom B
		Molib::Atom &atomA = seed.get_center_atom();
		Molib::Atom::Vec atoms;
		dbgmsg("atomA = " << atomA);
		
		// get seed atoms without atomA
		for (auto &patom : __ordered_atoms) {
			if (patom != &atomA) atoms.push_back(patom);
		}
		
		// get point that is at the center of the gpoints (0,0,0)
		Gpoints::Gpoint &pointA = const_cast<Gpoints::Gpoint&>(gpoints.get_center_point());
		vector<pair<Molib::Atom*, Gpoints::Gpoint*>> vertices;

		double coeff = 3.0;
		vector<vector<unsigned short int>> qmaxes;

		// iterate until you get a decent number of conformations
		do { 
			// generate product graph vertices
			vertices.clear();
			vertices.push_back(make_pair(&atomA, &pointA));
			const double tol = grid_spacing / coeff;
	
			for (auto &patomB : atoms) {
				for (auto &pointB : grid.get_neighbors_within_tolerance(pointA.crd(), atomA.crd().distance(patomB->crd()), tol)) {
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
						conn.set(i, j);
						conn.set(j, i);
					}
				}
			}
			
			// find all max cliques of size equal num. seed atoms
			Maxclique m(conn);
			qmaxes = m.mcq(__ordered_atoms.size());

			dbgmsg("iterative coeff = " << coeff << " qmaxes.size() = " << qmaxes.size());

			coeff -= 0.1;

		} while (qmaxes.size() < min_num_conf);
		
		// save conformations
		for (int c = 0; c < qmaxes.size(); ++c) {
			auto &qmax = qmaxes[c];
			// if max clique contains all seed atoms...
			dbgmsg("qmax[" << c << "].size = " << qmax.size());
			// go over vertices of each max clique
			Gpoints::PGpointVec points(qmax.size()); // init size of vec!
			for (auto &i : qmax) {
				Molib::Atom *patom = vertices[i].first;
				Gpoints::Gpoint *gpoint = vertices[i].second;
				Gpoints::IJK ijk = gpoint->ijk();
				__conf_map[ijk.i][ijk.j][ijk.k].push_back(c);
				points[atom_to_i[patom]] = gpoint;
			}
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
			<< seed.name() << " with coeff = " << coeff << " took " << Benchmark::seconds_from_start() 
			<< " wallclock seconds" << endl;
	}
	
};
