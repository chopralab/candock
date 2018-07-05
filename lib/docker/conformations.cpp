#include "candock/helper/inout.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/logger.hpp"
#include "candock/molib/grid.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/geometry/quaternion.hpp"
#include "candock/graph/mcqd.hpp"
#include "candock/helper/array2d.hpp"
#include "candock/docker/gpoints.hpp"
#include "candock/docker/conformations.hpp"
#include <iostream>
#include <exception>

using namespace std;

namespace candock{
namespace docker {

	ostream& operator<<(ostream& os, const Conformations &conformations) {
		const molib::Atom::Vec atoms = conformations.__seed.get_atoms();
		for (auto &conf : conformations.get_conformations()) {
			os << "MODEL" << endl;
			for (size_t i = 0; i < conf.size(); ++i) {
				Gpoints::Gpoint &gpoint0 = *conf[i];
				molib::Atom &atom = *atoms[i];
				atom.set_crd(gpoint0.crd());
				os << atom;
			}
			os << "ENDMDL" << endl;
		}
		return os;
	}

	Conformations::Conformations(const molib::Molecule &seed, const Gpoints &gpoints, 
		const double &conf_spin, const int num_univec) : __seed(seed) {

		try {
			const double conf_spin_in_radians = geometry::radians(conf_spin / 2);
			
			Benchmark bench;
	
			// get uniform sphere points, i.e., unit vectors
			geometry::Point::Vec unit_vectors = geometry::uniform_sphere(num_univec);
	
			const molib::Atom &center = seed.get_center_atom();
			
			// get another atom - neighbor of center
			const molib::Atom &another = center.first();
	
			// calculate vector between the two atoms
			const geometry::Point bondvec = another.crd() - center.crd();
			
			geometry::Point::Vec seed_crds = seed.get_crds();
	
			// translate points so that the central atom is at origin
			for (auto &crd : seed_crds)
				crd = crd - center.crd();
			
			// rotate seed on each unit_vector by increments of conf_spin_in_radians degrees (in radians)
			vector<geometry::Point::Vec> confs;
			for (auto &unitvec : unit_vectors) {
	
				geometry::Vector3 ortho = geometry::Coordinate::cross(bondvec, unitvec);
				const double rotangl = geometry::angle(unitvec, bondvec) / 2;
				dbgmsg("rotangl = " << geometry::degrees(rotangl) << " seed = " << seed.name());
				const geometry::Quaternion q0(ortho.norm()*sin(rotangl), cos(rotangl));
	
				// align seed crds along unit vector
				geometry::Point::Vec previous;
				for (auto &crd : seed_crds) {
					previous.push_back(q0.rotatedVector(crd)); 
				}

				confs.push_back(previous);
				
#ifndef NDEBUG
				stringstream ss;
				ss << "MODEL" << endl
					<< seed_crds 
					<< "ENDMDL" << endl
					<< "MODEL" << endl 
					<< previous 
					<< "ATOM      1   U  UNI     2    "  << unitvec.pdb() << endl 
					<< "ATOM      1   U  BON     3    "  << bondvec.pdb() << endl 
					<< "ENDMDL" << endl;
				Inout::output_file(ss.str(), "unit_" + seed.name() + ".pdb", ios_base::app); 
#endif				

				const geometry::Quaternion q(geometry::Vector3(unitvec)*sin(conf_spin_in_radians), cos(conf_spin_in_radians));
	
				for (double angle = 2 * conf_spin_in_radians; angle < M_PI; angle += conf_spin_in_radians) {
					geometry::Point::Vec rotated;
					for (auto &crd : previous) {	
						rotated.push_back(q.rotatedVector(crd)); 
					}
					confs.push_back(rotated);
					previous = rotated;
				}
			}
			// map rotamer conformations to gridpoints
			Gpoints::PGpointVec pgvec;
			for (auto &point : gpoints.get_gridpoints0()) {
				pgvec.push_back(const_cast<Gpoints::Gpoint*>(&point));
			 }
			molib::Grid<Gpoints::Gpoint> grid(pgvec);
			dbgmsg("after creating grid");
			dbgmsg("conformation size = " << confs.size());
	
			// pre-compute rmsd between ALL conformations
			__rmsd_sq.init(confs.size(), confs.size());
			for (size_t i = 0; i < confs.size(); ++i) {
				for (size_t j = i + 1; j < confs.size(); ++j) {
					__rmsd_sq.data[i][j] = __rmsd_sq.data[j][i] = geometry::compute_rmsd_sq(confs[i], confs[j]);
				}
			}
			dbgmsg("after rmsd calculation");
	
			// go over all conformations & get the closest gripoint
			for (size_t i = 0; i < confs.size(); ++i) {
	
				auto &conf = confs[i];
				Gpoints::PGpointVec points; // init size of vec!
				dbgmsg("conformation size = " << conf.size());
				for (auto &crd : conf) {
					Gpoints::Gpoint *closest = nullptr;
					double max_dist = HUGE_VAL;
					for (auto &pgpoint : grid.get_neighbors(crd, 1.0)) {
						if (pgpoint->crd().distance(crd) < max_dist) {
							closest = pgpoint;
							max_dist = pgpoint->crd().distance(crd);
						}
					}
					if (!closest) throw Error("die : fragment is outside of grid");
					Gpoints::IJK ijk = closest->ijk();
					__conf_map[ijk.i][ijk.j][ijk.k].push_back(i);
					points.push_back(closest);
	
				}
				// save conformations
				__conf_vec.push_back(points);
			}

                        log_benchmark << "time to find " << __conf_vec.size() << " conformations of seed " 
                                      << seed.name() << " took " << bench.seconds_from_start() 
                                      << " wallclock seconds" << "\n";
	
	
		} catch(...) {
			dbgmsg("FAILURE: constructor of Conformations failed ... cleaning resources...");
			throw;
		}
	}

	
}
}
