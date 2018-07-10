#include "candock/helper/inout.hpp"
#include "candock/molib/grid.hpp"
#include "candock/score/score.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/debug.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/helper/array1d.hpp"
#include "candock/docker/gpoints.hpp"
#include "candock/docker/conformations.hpp"
#include "candock/docker/dock.hpp"
#include <iostream>
#include <exception>
using namespace std;

namespace candock{
namespace docker {

        double Dock::DockedConf::compute_rmsd_sq(const Dock::DockedConf &other) const {
        
                const Gpoints::PGpointVec &points1 = this->get_conf0();
                const Gpoints::PGpointVec &points2 = other.get_conf0();

                double rmsd = 0.0;

                for ( size_t i = 0; i < points1.size(); ++i ) {
                        rmsd += points1[i]->crd().distance_sq(points2[i]->crd());
                }

                return rmsd / points1.size();
        }

	void Dock::run() {
		DockedConf::Vec docked = __dock();
		DockedConf::Vec clustered = __cluster(docked);
		docked.clear(); // clear memory
		__set_docked(clustered);
	}

	Dock::DockedConf::Vec Dock::__dock() {
		Benchmark bench;

		DockedConf::Vec accepted;
		
		auto &conformations = __conformations.get_conformations();
		Array1d<bool> rejected(conformations.size());

		molib::Atom::Vec seed_atoms = __seed.get_atoms();

		// go over all cavity points
		for (auto &kv : __gpoints.get_gridpoints()) {
            const int bsite_id = kv.first;
            auto &gmap = __gpoints.get_gmap(bsite_id);
			for (auto &cavpoint : kv.second) {
				DockedConf::Vec accepted_tmp;
				// reset map of rejected conformations to zero
				rejected.reset();
				for (size_t c = 0; c < conformations.size(); ++c) {
					Gpoints::PGpointVec &conf = conformations[c];
					// test if c-th conformation clashes with receptor: if yes, reject it
					if (!rejected.data[c]) {
						double energy_sum(0);
						bool reje = false;
						// go over coordinates of the c-th conformation
						dbgmsg("testing conformation " << c);

						for (size_t i = 0; i < conf.size(); ++i) {
							docker::Gpoints::Gpoint &gpoint0 = *conf[i];
							docker::Gpoints::IJK confijk = cavpoint.ijk() + gpoint0.ijk();

							dbgmsg("cavpoint.ijk() = " << cavpoint.ijk());
							dbgmsg("gpoint.ijk() = " << gpoint0.ijk());
							dbgmsg("confijk = " << confijk);
							dbgmsg("gmap.szi = " << gmap.szi << " gmap.szj = " << gmap.szj
								<< " gmap.szk = " << gmap.szk);


							if (confijk.i < 0 || confijk.j < 0 || confijk.k < 0 ||
								(size_t)confijk.i >= gmap.szi || (size_t)confijk.j >= gmap.szj || (size_t)confijk.k >= gmap.szk
								|| gmap.data[confijk.i][confijk.j][confijk.k] == nullptr) {
								// mark as rejected all conformations that have this point
								for (auto &r : __conformations.get_confs_at(gpoint0.ijk())) {
									rejected.data[r] = true;
									dbgmsg("rejected conformation " << r);
								}
								reje = true;
								break;
							}
							docker::Gpoints::Gpoint *pgpoint = gmap.data[confijk.i][confijk.j][confijk.k];
							molib::Atom &atom = *seed_atoms[i];
							energy_sum += pgpoint->energy(atom.idatm_type());
							dbgmsg("energy of pgpoint at crd = " << pgpoint->crd() 
								<< " is " << pgpoint->energy(atom.idatm_type()) 
								<< " for idatm_type = " << atom.idatm_type());
							dbgmsg("energy calculated at crd = " << pgpoint->crd() 
								<< " is " << __score.non_bonded_energy(__gridrec, 
								molib::Atom::Vec{&atom}, geometry::Point::Vec{pgpoint->crd()}));
						}
						// if no clashes were found ...
						if (!reje) {
							accepted_tmp.push_back(Dock::DockedConf(cavpoint, conf, energy_sum, c, bsite_id));
#ifndef NDEBUG
							molib::Atom::Vec seed_atoms = __seed.get_atoms();
							const Gpoints::PGpointVec &points = accepted_tmp.back().get_conf0();

							for (size_t i = 0; i < points.size(); ++i) {
								
								molib::Atom &atom = *seed_atoms[i];
								Docker::Gpoints::Gpoint &gpoint0 = *points[i];
								Docker::Gpoints::IJK confijk = accepted_tmp.back().get_cavpoint().ijk() + gpoint0.ijk();
								Docker::Gpoints::Gpoint *pgpoint = gmap.data[confijk.i][confijk.j][confijk.k];
								
								atom.set_crd(pgpoint->crd());
								dbgmsg("crd = " << atom.crd());
								
							}
							dbgmsg("ACCEPTED conformation energy = " << accepted_tmp.back().get_energy() 
								<< " calculated energy = " << __score.non_bonded_energy(__gridrec, __seed));
#endif
						}
					}
				}
				__cluster_fast(accepted_tmp, accepted);
			}
		}
		log_benchmark << "Fragment docking for seed " << __seed.name() << " took " << bench.seconds_from_start() << " seconds"
			<< " number of accepted confomations is " << accepted.size() << "\n";

		return accepted;
	}

	void Dock::__cluster_fast(const DockedConf::Vec &conformations, DockedConf::Vec &reps) {

		set<const Dock::DockedConf*, Dock::DockedConf::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&conf);

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const Dock::DockedConf &lowest_point = **confs.begin();
			reps.push_back(lowest_point);
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy conformation
			for (auto it = confs.begin(); it != confs.end(); ++it) {
				if (__conformations.get_rmsd_sq(lowest_point.get_i(), (*it)->get_i()) < __rmsd_tol_sq) {
					confs.erase(it);
				}
			}
		}
	}

	Dock::DockedConf::Vec Dock::__cluster(const DockedConf::Vec &conformations) {

		Benchmark bench;
		DockedConf::Vec reps;

		set<const Dock::DockedConf*, Dock::DockedConf::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&conf);
	
		molib::Grid<const Dock::DockedConf> cgrid(confs); // grid of docked conformations

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const Dock::DockedConf &lowest_point = **confs.begin();
			reps.push_back(lowest_point);
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy yconformation
			for (auto &pconf : cgrid.get_neighbors(lowest_point.crd(), __rmsd_tol)) {
				if (pconf->compute_rmsd_sq(lowest_point) < __rmsd_tol_sq) {
					confs.erase(pconf);
				}
			}
		}
		log_benchmark << "Clustering accepted conformations for seed " << __seed.name() << " took " 
			<< bench.seconds_from_start() << " seconds" << "\n";
		return reps;
	}

	void Dock::__set_docked(const DockedConf::Vec &confs) {

		Benchmark bench;

		__docked.set_name(__seed.name()); // molecules(!) name is seed_id
		
		molib::Atom::Vec seed_atoms = __seed.get_atoms();

		// go over all accepted conformations
		for (auto &conf : confs) {
			auto &gmap = __gpoints.get_gmap(conf.get_bsite_id());

			dbgmsg(" conformation size = " << conf.get_conf0().size() 
				<< " energy = " << conf.get_energy());
			// correct the seed's new coordinates and ...
			const Gpoints::PGpointVec &points = conf.get_conf0();
			for (size_t i = 0; i < points.size(); ++i) {

				molib::Atom &atom = *seed_atoms[i];
				docker::Gpoints::Gpoint &gpoint0 = *points[i];
	
				docker::Gpoints::IJK confijk = conf.get_cavpoint().ijk() + gpoint0.ijk();
				docker::Gpoints::Gpoint *pgpoint = gmap.data[confijk.i][confijk.j][confijk.k];

				atom.set_crd(pgpoint->crd());
			}
			// save the conformation
			__docked.add(new molib::Molecule(__seed)).set_name(std::to_string(conf.get_energy())); // molecule name is energy
			dbgmsg("conformation energy = " << conf.get_energy() 
				<< " calculated energy = " << __score.non_bonded_energy(__gridrec, __seed));

		}

		dbgmsg("Conversion of conformations to mols took " 
			<< bench.seconds_from_start() << " seconds");
	}

}
}
