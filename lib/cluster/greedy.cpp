#include "candock/molib/grid.hpp"
#include "candock/molib/molecules.hpp"
#include "candock/score/score.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/helper/debug.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/cluster/greedy.hpp"
#include <iostream>
#include <exception>
#include "candock/helper/logger.hpp"

using namespace candock::molib;

namespace candock{
namespace Cluster {

	/**
	 * Cluter points
	 */
	geometry::Point::Vec Cluster::greedy(const geometry::Point::Vec &initial, const double clus_rad) {

		Benchmark bench;

		geometry::Point::ConstSet confs;
		for (auto &point : initial)
			confs.insert(&point);

		dbgmsg(confs);
	
		Grid<const geometry::Point> grid(confs); // grid of points

		geometry::Point::Vec reps;

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const geometry::Point &lowest_point = **confs.begin();
			reps.push_back(lowest_point);
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy yconformation
			for (auto &pconf : grid.get_neighbors(lowest_point.crd(), clus_rad)) {
				confs.erase(pconf);
			}
		}
		log_benchmark << "Clustering " << initial.size() << " accepted conformations resulted in "
			<< reps.size() << " clusters took " << bench.seconds_from_start() 
			<< " seconds" << "\n";
		return reps;
	}

	molib::Molecules Cluster::greedy(const molib::Molecules &initial, const score::Score &score, molib::Atom::Grid &gridrec, const double clus_rad) {

		Benchmark bench;

		LinkedConf<Molecule>::UPVec conformations;
		for (auto &molecule : initial)
			conformations.push_back(unique_ptr<LinkedConf<Molecule>>(new LinkedConf<Molecule>(molecule, molecule.compute_geometric_center(), score.non_bonded_energy(gridrec, molecule))));

		set<const LinkedConf<Molecule>*, LinkedConf<Molecule>::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&*conf);

		dbgmsg(confs);
	
		Grid<const LinkedConf<Molecule>> grid(confs); // grid of docked conformations

		molib::Molecules reps;

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const LinkedConf<Molecule> &lowest_point = **confs.begin();
			reps.add(new molib::Molecule(lowest_point.get_molecule()));
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy yconformation
			for (auto &pconf : grid.get_neighbors(lowest_point.crd(), clus_rad)) {
				if (pconf->get_molecule().compute_rmsd_ord(lowest_point.get_molecule()) < clus_rad) {
					confs.erase(pconf);
				}
			}
		}
		log_benchmark << "Clustering " << initial.size() << " accepted conformations resulted in "
			<< reps.size() << " clusters took " << bench.seconds_from_start() 
			<< " seconds" << "\n";
		return reps;
	}

	linker::Partial::Vec Cluster::greedy(const linker::Partial::Vec &initial, const molib::Atom::Grid&, const double clus_rad) {

		Benchmark bench;

		LinkedConf<linker::Partial>::UPVec conformations;
		for (auto &molecule : initial) {
			conformations.push_back(unique_ptr<LinkedConf<linker::Partial>>(new LinkedConf<linker::Partial>(const_cast<linker::Partial&>(molecule), molecule.compute_geometric_center(), molecule.get_energy())));
		}

		set<const LinkedConf<linker::Partial>*, LinkedConf<linker::Partial>::by_energy> confs;
		for (auto &conf : conformations) confs.insert(&*conf);

		dbgmsg("number of conformations to cluster = " << confs.size() << endl << confs);
	
		Grid<const LinkedConf<linker::Partial>> grid(confs); // grid of docked conformations

		linker::Partial::Vec reps;

		while (!confs.empty()) {
			// accept lowest energy conformation as representative
			const LinkedConf<linker::Partial> &lowest_point = **confs.begin();
			reps.push_back(lowest_point.get_molecule());
			confs.erase(confs.begin());
			// delete all conformations within RMSD tolerance of this lowest energy conformation
			for (auto &pconf : grid.get_neighbors(lowest_point.crd(), clus_rad)) {
				try {
					dbgmsg("found neighbor of conformation");
					if (pconf->get_molecule().compute_rmsd_ord(lowest_point.get_molecule()) < clus_rad) {
						dbgmsg("rmsd between two conformations is " 
							<< pconf->get_molecule().compute_rmsd_ord(lowest_point.get_molecule()));
						dbgmsg("conformation1 " << lowest_point.get_molecule());
						dbgmsg("conformation2 " << pconf->get_molecule());
						confs.erase(pconf);
					}
				} catch(const Error &) {
					// if calculation of rmsd didn't succeed don't do nothing
					// consequently all conformations will be representative
				}
			}
		}
		log_benchmark << "Clustering " << initial.size() << " accepted conformations resulted in "
			<< reps.size() << " clusters took " << bench.seconds_from_start() 
			<< " seconds" << "\n";
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

	ostream& operator<<(ostream& os, const set<const Cluster::LinkedConf<linker::Partial>*, Cluster::LinkedConf<linker::Partial>::by_energy> &confs)	{
		for (auto &pconf : confs) {
			os << "linked conformation : " << " crd = " << pconf->crd() 
				<< " energy = "	<< setprecision(2) 
				<< fixed << pconf->get_energy() << endl;
		}
		return os;
	}

}
}
