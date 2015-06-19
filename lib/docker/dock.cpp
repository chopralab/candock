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

	//~ Molib::Molecules Dock::run() {
	AcceptedConformations Dock::run() {
		Benchmark::reset();

		//~ Molib::Molecules docked;
		AcceptedConformations accepted;
		
		auto cavity_points = __gpoints.get_gridpoints_as_vec();
		auto &gmap = __gpoints.get_gmap();
		auto &conformations = __conformations.get_conformations();
		auto &confmap = __conformations.get_confmap();
		Array1d<bool> rejected(conformations.size());

		int acc_i = 0;
		
		// go over all cavity points
		for (auto &cavpoint : cavity_points) {
			// reset map of rejected conformations to 0
			rejected.reset();
			for (int c = 0; c < conformations.size(); ++c) {
				// test if c-th conformation clashes with receptor: if yes, reject it
				if (!rejected.data[c]) {
					// go over coordinates of the c-th conformation
					dbgmsg("testing conformation " << c);
					for (auto &pair : conformations[c]) {
						Docker::Gpoint &gpoint = *pair.second;
						Docker::IJK confijk = cavpoint.ijk() + gpoint.ijk();
						dbgmsg("cavpoint.ijk() = " << cavpoint.ijk());
						dbgmsg("gpoint.ijk() = " << gpoint.ijk());
						dbgmsg("confijk = " << confijk);
						dbgmsg("gmap.szi = " << gmap.szi << " gmap.szj = " << gmap.szj
							<< " gmap.szk = " << gmap.szk);
						if (confijk.i < 0 || confijk.j < 0 || confijk.k < 0 ||
							confijk.i >= gmap.szi || confijk.j >= gmap.szj || confijk.k >= gmap.szk
							|| gmap.data[confijk.i][confijk.j][confijk.k] == nullptr) {
							// mark as rejected all conformations that have this point
							for (auto &cr : confmap[gpoint.ijk().i][gpoint.ijk().j][gpoint.ijk().k]) {
								rejected.data[cr] = true;
								dbgmsg("rejected conformation " << cr);
							}
							break;
						}
					}
					// if no clashes were found ...
					if (!rejected.data[c]) {
						++acc_i;
						accepted.push_back({&cavpoint, &conformations[c]});
						dbgmsg("conformation accepted");
					}
				}
			}
		}
		cout << "Fragment docking took " << Benchmark::seconds_from_start() << " seconds"
			<< " number of accepted confomations is " << acc_i << endl;

		//~ return docked;
		return accepted;
	}

	Molib::Molecules Dock::cluster(const AcceptedConformations &accepted) {

		Benchmark::reset();

		Molib::Molecules docked;
		auto &gmap = __gpoints.get_gmap();

		// go over all accepted conformations
		for (auto &ac : accepted) {
			Gpoint &cavpoint = *ac.first;
			OneConformation &conf = *ac.second;
			double energy_sum = 0;
			// score this conformation
			for (auto &pair : conf) {
				Molib::Atom &atom = *pair.first;
				Docker::Gpoint &gpoint0 = *pair.second;
				Docker::IJK confijk = cavpoint.ijk() + gpoint0.ijk();
				Docker::Gpoint &gpoint = *gmap.data[confijk.i][confijk.j][confijk.k];
				energy_sum += gpoint.energy(atom.idatm_type());
			}
			
			//~ if (energy_sum > 100000) {
				//~ cout << " conf.size() = " << conf.size() << " energy_sum = " << energy_sum << endl;			
				//~ // correct the seed's new coordinates and ...
				//~ cout << "MODEL" << endl;
				//~ for (auto &pair : conf) {
					//~ Molib::Atom &atom = *(pair.first);
					//~ Docker::Gpoint &gpoint0 = *pair.second;
					//~ atom.set_crd(cavpoint.crd() + gpoint0.crd());
					//~ cout << atom;
				//~ }
				//~ cout << "ENDMDL" << endl;
				//~ // save the c-th conformation
				//~ docked.add(new Molib::Molecule(__seed));
			//~ }
		}

		//~ Benchmark::reset();
		//~ Geom3D::PointVec v;
		//~ for (int i = 0; i < accepted.size(); ++i) {
			//~ Gpoint &cavpoint = *accepted[i].first;
			//~ OneConformation &conf = *accepted[i].second;
			//~ for (int j = i + 1; j < accepted.size(); ++j) {
				//~ Gpoint &cavpoint2 = *accepted[j].first;
				//~ OneConformation &conf2 = *accepted[j].second;
				//~ if (cavpoint.crd().distance(cavpoint2.crd()) < 0.5) {
					//~ v.push_back(cavpoint.crd());
				//~ }
			//~ }
		//~ }
//~ 
		//~ cout << "v.size() = " << v.size() << endl;
		
		cout << "Clustering accepted conformations took " << Benchmark::seconds_from_start() 
			<< " seconds" << endl;

		return docked;
	}
};
