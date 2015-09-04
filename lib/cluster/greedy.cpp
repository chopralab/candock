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

		LinkedConf<Molecule>::UPVec conformations;
		for (auto &molecule : initial)
			conformations.push_back(unique_ptr<LinkedConf<Molecule>>(new LinkedConf<Molecule>(molecule, molecule.compute_geometric_center(), score.non_bonded_energy(gridrec, molecule))));

		set<const LinkedConf<Molecule>*, LinkedConf<Molecule>::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&*conf);

		dbgmsg(confs);
	
		Grid<const LinkedConf<Molecule>> grid(confs); // grid of docked conformations

		Molib::Molecules reps;

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const LinkedConf<Molecule> &lowest_point = **confs.begin();
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

	Linker::Partial::Vec Cluster::greedy(const Linker::Partial::Vec &initial, Molib::Atom::Grid &gridrec, const double clus_rad) {

		Benchmark::reset();

		LinkedConf<Linker::Partial>::UPVec conformations;
		for (auto &molecule : initial) {
			conformations.push_back(unique_ptr<LinkedConf<Linker::Partial>>(new LinkedConf<Linker::Partial>(const_cast<Linker::Partial&>(molecule), molecule.compute_geometric_center(), molecule.get_energy())));
		}

		set<const LinkedConf<Linker::Partial>*, LinkedConf<Linker::Partial>::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&*conf);

		dbgmsg("number of conformations to cluster = " << confs.size() << endl << confs);
	
		Grid<const LinkedConf<Linker::Partial>> grid(confs); // grid of docked conformations

		Linker::Partial::Vec reps;

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const LinkedConf<Linker::Partial> &lowest_point = **confs.begin();
			reps.push_back(lowest_point.get_molecule());
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

	ostream& operator<<(ostream& os, const set<const Cluster::LinkedConf<Molecule>*, Cluster::LinkedConf<Molecule>::by_energy> &confs)	{
		for (auto &pconf : confs) {
			os << "docked conformation : name = " << pconf->get_molecule().name() 
				<< " crd = " << pconf->crd() << " energy = "	<< setprecision(2) 
				<< fixed << pconf->get_energy() << endl;
		}
		return os;
	}

	ostream& operator<<(ostream& os, const set<const Cluster::LinkedConf<Linker::Partial>*, Cluster::LinkedConf<Linker::Partial>::by_energy> &confs)	{
		for (auto &pconf : confs) {
			os << "linked conformation : " << " crd = " << pconf->crd() 
				<< " energy = "	<< setprecision(2) 
				<< fixed << pconf->get_energy() << endl;
		}
		return os;
	}

};
