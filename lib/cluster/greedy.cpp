#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "score/score.hpp"
#include "helper/benchmark.hpp"
#include "helper/debug.hpp"
#include "geom3d/geom3d.hpp"
#include "greedy.hpp"
#include <iostream>
#include <exception>

namespace Molib {

	Molib::Molecules Cluster::greedy(const Molib::Molecules &initial, const Molib::Score &score, Molib::Atom::Grid &gridrec, const double clus_rad) {

		Benchmark::reset();

		LinkedConf::UPVec conformations;
		for (auto &molecule : initial)
			conformations.push_back(unique_ptr<LinkedConf>(new LinkedConf(molecule, molecule.compute_geometric_center(), score.non_bonded_energy(gridrec, molecule))));

		set<const LinkedConf*, LinkedConf::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&*conf);

		dbgmsg(confs);
	
		Grid<const LinkedConf> grid(confs); // grid of docked conformations

		Molib::Molecules reps;

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const LinkedConf &lowest_point = **confs.begin();
			reps.add(new Molib::Molecule(lowest_point.get_molecule()));
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy yconformation
			for (auto &pconf : grid.get_neighbors(lowest_point.crd(), clus_rad)) {
				if (pconf->get_molecule().compute_rmsd_ord(lowest_point.get_molecule()) < clus_rad) {
					confs.erase(pconf);
				}
			}
		}
		cout << "Clustering " << initial.size() << " accepted conformations resulted in "
			<< reps.size() << " clusters took " << Benchmark::seconds_from_start() 
			<< " seconds" << endl;
		return reps;
	}

	ostream& operator<<(ostream& os, const set<const Cluster::LinkedConf*, Cluster::LinkedConf::by_energy> &confs)	{
		for (auto &pconf : confs) {
			os << "docked conformation : name = " << pconf->get_molecule().name() 
				<< " crd = " << pconf->crd() << " energy = "	<< setprecision(2) 
				<< fixed << pconf->get_energy() << endl;
		}
		return os;
	}

};
