#include "helper/inout.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "score/score.hpp"
#include "helper/benchmark.hpp"
#include "helper/debug.hpp"
#include "geom3d/geom3d.hpp"
#include "helper/array1d.hpp"
#include "gpoints.hpp"
#include "conformations.hpp"
#include "dock.hpp"
#include <iostream>
#include <exception>

namespace Docker {

	double Dock::DockedConf::compute_rmsd(const Dock::DockedConf &other) const {
		
		double sum_sq(0);
		for (int i = 0; i < this->__conf0.get_points().size(); ++i) {
			Geom3D::Point crda = this->__cavpoint.crd() + this->__conf0.get_point(i).crd();
			Geom3D::Point crdb = other.__cavpoint.crd() + other.__conf0.get_point(i).crd();
			sum_sq += crda.distance_sq(crdb);
		}
		return sqrt(sum_sq / this->__conf0.get_points().size());
	}

	Molib::Molecules Dock::run() {
		DockedConfVec docked = __dock();
		DockedConfVec clustered = __cluster(docked);
		return __convert_to_mols(clustered);
	}

	Dock::DockedConfVec Dock::__dock() {
		Benchmark::reset();

		DockedConfVec accepted;
		
		auto &gmap = __gpoints.get_gmap();
		auto &conformations = __conformations.get_conformations();
		auto &confmap = __conformations.get_confmap();
		Array1d<bool> rejected(conformations.size());

		// go over all cavity points
		for (auto &kv : __gpoints.get_gridpoints()) {
			for (auto &cavpoint : kv.second) {
				DockedConfVec accepted_tmp;
				// reset map of rejected conformations to zero
				rejected.reset();
				for (int c = 0; c < conformations.size(); ++c) {
					Conformations::Conf &conf = conformations[c];
					// test if c-th conformation clashes with receptor: if yes, reject it
					if (!rejected.data[c]) {
						double energy_sum(0);
						// go over coordinates of the c-th conformation
						dbgmsg("testing conformation " << c);
						for (int i = 0; i < conf.get_points().size(); ++i) {
							Docker::Gpoints::Gpoint &gpoint0 = conf.get_point(i);
							Docker::Gpoints::IJK confijk = cavpoint.ijk() + gpoint0.ijk();

							dbgmsg("cavpoint.ijk() = " << cavpoint.ijk());
							dbgmsg("gpoint.ijk() = " << gpoint0.ijk());
							dbgmsg("confijk = " << confijk);
							dbgmsg("gmap.szi = " << gmap.szi << " gmap.szj = " << gmap.szj
								<< " gmap.szk = " << gmap.szk);
							if (confijk.i < 0 || confijk.j < 0 || confijk.k < 0 ||
								confijk.i >= gmap.szi || confijk.j >= gmap.szj || confijk.k >= gmap.szk
								|| gmap.data[confijk.i][confijk.j][confijk.k] == nullptr) {
								// mark as rejected all conformations that have this point
								for (auto &r : __conformations.get_confs_at(gpoint0.ijk())) {
									rejected.data[r] = true;
									dbgmsg("rejected conformation " << r);
								}
								goto ESCAPE;
							}
							Docker::Gpoints::Gpoint *pgpoint = gmap.data[confijk.i][confijk.j][confijk.k];
							Molib::Atom &atom = conf.get_atom(i);
							energy_sum += pgpoint->energy(atom.idatm_type());
						}
						// if no clashes were found ...
						accepted_tmp.push_back(Dock::DockedConf(cavpoint, conf, energy_sum, c));
						
						ESCAPE:
						;
					}
				}
				__cluster_fast(accepted_tmp, accepted);
			}
		}
		cout << "Fragment docking took " << Benchmark::seconds_from_start() << " seconds"
			<< " number of accepted confomations is " << accepted.size() << endl;

		return accepted;
	}

	void Dock::__cluster_fast(const DockedConfVec &conformations, DockedConfVec &reps) {

		set<const Dock::DockedConf*, Dock::DockedConf::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&conf);

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const Dock::DockedConf &lowest_point = **confs.begin();
			reps.push_back(lowest_point);
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy conformation
			for (auto it = confs.begin(); it != confs.end(); ++it) {
				if (__conformations.get_rmsd(lowest_point.get_i(), (*it)->get_i()) < __rmsd_tol) {
					confs.erase(it);
				}
			}
		}
	}

	Dock::DockedConfVec Dock::__cluster(const DockedConfVec &conformations) {

		Benchmark::reset();
		DockedConfVec reps;

		set<const Dock::DockedConf*, Dock::DockedConf::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&conf);
	
		Grid<const Dock::DockedConf> cgrid(confs); // grid of docked conformations

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const Dock::DockedConf &lowest_point = **confs.begin();
			reps.push_back(lowest_point);
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy yconformation
			for (auto &pconf : cgrid.get_neighbors(lowest_point.crd(), __rmsd_tol)) {
				if (pconf->compute_rmsd(lowest_point) < __rmsd_tol) {
					confs.erase(pconf);
				}
			}
		}
		cout << "Clustering accepted conformations took " 
			<< Benchmark::seconds_from_start() << " seconds" << endl;
		return reps;
	}

	Molib::Molecules Dock::__convert_to_mols(const DockedConfVec &confs) {

		Benchmark::reset();
		
		Molib::Molecules docked;
		docked.set_name(__seed.name()); // molecules name is seed_id
		// go over all accepted conformations
		for (auto &conf : confs) {
			dbgmsg(" conformation size = " << conf.get_conf0().get_points().size() 
				<< " energy = " << conf.get_energy());
			// correct the seed's new coordinates and ...
			for (int i = 0; i < conf.get_conf0().get_points().size(); ++i) {
				Molib::Atom &atom = conf.get_conf0().get_atom(i);
				Docker::Gpoints::Gpoint &gpoint0 = conf.get_conf0().get_point(i);
				atom.set_crd(conf.get_cavpoint().crd() + gpoint0.crd());
			}
			//~ __seed.set_name(help::to_string(conf.get_energy())); // molecule name is energy
			// save the conformation
			//~ docked.add(new Molib::Molecule(__seed));
			docked.add(new Molib::Molecule(__seed)).set_name(help::to_string(conf.get_energy())); // molecule name is energy

		}
		dbgmsg("Conversion of conformations to mols took " 
			<< Benchmark::seconds_from_start() << " seconds");
		return docked;
	}

};
