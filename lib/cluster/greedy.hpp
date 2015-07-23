#ifndef GREEDY_CLUSTER_H
#define GREEDY_CLUSTER_H

#include "pdbreader/molecule.hpp"
#include "geom3d/geom3d.hpp"

namespace Molib {
	class Score;
	class Molecule;
	class Molecules;
	
	class Cluster {
		class LinkedConf {
		public:
			typedef vector<unique_ptr<LinkedConf>> UPVec;
		private:
			Molib::Molecule &__molecule;
			Geom3D::Point __crd;
			double __energy;
		public:
			LinkedConf(Molib::Molecule &molecule, Geom3D::Point crd, double energy) 
				: __molecule(molecule), __crd(crd), __energy(energy) {}
			Geom3D::Point& crd() { return __crd; }
			const Geom3D::Point& crd() const { return __crd; }
			void distance(const double) const {} // dummy
			Molib::Molecule &get_molecule() const { return __molecule; }
			double get_energy() const { return __energy; }
			struct by_energy {
				bool operator() (const LinkedConf *lhs, const LinkedConf *rhs) const {
					return lhs->__energy < rhs->__energy;
				}
			};
		};
		friend ostream& operator<<(ostream& os, const set<const LinkedConf*, LinkedConf::by_energy> &confs);
	public:
		static Molib::Molecules greedy(const Molib::Molecules &initial, const Molib::Score &score,
			const double clus_rad);
	};
		
};

#endif
