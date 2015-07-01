#ifndef DOCK_H
#define DOCK_H

#include "gpoints.hpp"
#include "conformations.hpp"

namespace Molib {
	class Molecule;
};

namespace Docker {
	
	class Dock {

		class DockedConf {
			Gpoints::Gpoint &__cavpoint;
			Conformations::Conf &__conf0;
			double __energy;
			int __i;
		public:
			DockedConf(Gpoints::Gpoint &cavpoint, Conformations::Conf &conf0, double energy, int i) 
				: __cavpoint(cavpoint), __conf0(conf0), __energy(energy), __i(i) {}
			Geom3D::Point& crd() { return __cavpoint.crd(); }
			const Geom3D::Point& crd() const { return __cavpoint.crd(); }
			void distance(const double) const {} // dummy
			Gpoints::Gpoint &get_cavpoint() const { return __cavpoint; }
			Conformations::Conf &get_conf0() const { return __conf0; }
			double get_energy() const { return __energy; }
			int get_i() const { return __i; }
			struct by_energy {
				bool operator() (const DockedConf *lhs, const DockedConf *rhs) const {
					return lhs->__energy < rhs->__energy;
				}
			};			
			double compute_rmsd(const DockedConf &other) const;
		};

		typedef vector<DockedConf> DockedConfVec;
		
		Gpoints &__gpoints;
		Conformations &__conformations;
		Molib::Molecule &__seed;
		double __rmsd_tol;

		DockedConfVec __dock();
		DockedConfVec __cluster(const DockedConfVec &confs);
		void __cluster_fast(const DockedConfVec &conformations, DockedConfVec &reps);
		Molib::Molecules __convert_to_mols(const DockedConfVec &confs);
		
	public:
		Dock(Gpoints &gpoints, Conformations &conformations,
			Molib::Molecule &seed, const double rmsd_tol=2.0) : __gpoints(gpoints),
			__conformations(conformations), __seed(seed), __rmsd_tol(rmsd_tol) {}
		Molib::Molecules run();
	};
};

#endif
