#ifndef DOCK_H
#define DOCK_H

#include "gpoints.hpp"
#include "conformations.hpp"
#include "pdbreader/molecules.hpp"

namespace Molib {
	class Molecule;
};

namespace Docker {
	
	class Dock {

		class DockedConf {
		public:
			typedef vector<DockedConf> Vec;

		private:
			const Gpoints::Gpoint &__cavpoint;
			Gpoints::PGpointVec &__conf0;
			double __energy;
			int __i;
		public:
			DockedConf(const Gpoints::Gpoint &cavpoint, Gpoints::PGpointVec &conf0, double energy, int i) 
				: __cavpoint(cavpoint), __conf0(conf0), __energy(energy), __i(i) {}
			//Geom3D::Point& crd() { return __cavpoint.crd(); }
			const Geom3D::Point& crd() const { return __cavpoint.crd(); }
			void distance(const double) const {} // dummy
			const Gpoints::Gpoint &get_cavpoint() const { return __cavpoint; }
			const Gpoints::PGpointVec &get_conf0() const { return __conf0; }
			double get_energy() const { return __energy; }
			int get_i() const { return __i; }
			struct by_energy {
				bool operator() (const DockedConf *lhs, const DockedConf *rhs) const {
					return lhs->__energy < rhs->__energy;
				}
			};			
			double compute_rmsd(const DockedConf &other) const;
		};

		const Gpoints &__gpoints;
		Conformations &__conformations;
		const Molib::Molecule &__seed;
		double __rmsd_tol;

		Molib::Molecules __docked;

#ifndef NDEBUG
		const Molib::Score &__score;
		const Molib::Atom::Grid &__gridrec;
#endif
		DockedConf::Vec __dock();
		DockedConf::Vec __cluster(const DockedConf::Vec &confs);
		void __cluster_fast(const DockedConf::Vec &conformations, DockedConf::Vec &reps);
		void __set_docked(const DockedConf::Vec &confs);
		
	public:
#ifndef NDEBUG
		Dock(const Gpoints &gpoints, Conformations &conformations, const Molib::Molecule &seed, const Molib::Score &score, 
			 const Molib::Atom::Grid &gridrec, const double rmsd_tol=2.0) : __gpoints(gpoints), __conformations(conformations), 
			__seed(seed), __score(score), __gridrec(gridrec), __rmsd_tol(rmsd_tol) {}
#else
		Dock(const Gpoints &gpoints, Conformations &conformations, Molib::Molecule &seed, const double rmsd_tol=2.0) 
			: __gpoints(gpoints), __conformations(conformations), __seed(seed), __rmsd_tol(rmsd_tol) {}
#endif
		void run();
		Molib::Molecules& get_docked() { return __docked; };
	};
};

#endif
