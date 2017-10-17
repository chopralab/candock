#ifndef GREEDY_CLUSTER_H
#define GREEDY_CLUSTER_H

#include "geom3d/geom3d.hpp"
#include "linker/partial.hpp"
#include "score/score.hpp"

namespace Score {
        class Score;
}

namespace Molib {
        class Molecule;
        class Molecules;
}

namespace Cluster {
	class Cluster {
		template<typename T>
		class LinkedConf {
		public:
			typedef vector<unique_ptr<LinkedConf>> UPVec;
		private:
			T &__molecule;
			Geom3D::Point __crd;
			double __energy;
		public:
			LinkedConf(T &molecule, Geom3D::Point crd, double energy) 
				: __molecule(molecule), __crd(crd), __energy(energy) {}
			Geom3D::Point& crd() { return __crd; }
			const Geom3D::Point& crd() const { return __crd; }
			void distance(const double) const {} // dummy
			T &get_molecule() const { return __molecule; }
			double get_energy() const { return __energy; }
			struct by_energy {
				bool operator() (const LinkedConf *lhs, const LinkedConf *rhs) const {
					return lhs->__energy < rhs->__energy;
				}
			};
		};
		friend ostream& operator<<(ostream& os, const set<const LinkedConf<Molib::Molecule>*, LinkedConf<Molib::Molecule>::by_energy> &confs);
		friend ostream& operator<<(ostream& os, const set<const LinkedConf<Linker::Partial>*, LinkedConf<Linker::Partial>::by_energy> &confs);
	public:
		static Geom3D::Point::Vec greedy(const Geom3D::Point::Vec &initial, const double clus_rad);
		static Molib::Molecules greedy(const Molib::Molecules &initial, const Score::Score &score, Molib::Atom::Grid &gridrec, const double clus_rad);
		static Linker::Partial::Vec greedy(const Linker::Partial::Vec &initial, const Molib::Atom::Grid &gridrec, const double clus_rad);
	};
		
};

#endif
