#ifndef GREEDY_CLUSTER_H
#define GREEDY_CLUSTER_H

#include "candock/geometry/geometry.hpp"
#include "candock/linker/partial.hpp"
#include "candock/score/score.hpp"

namespace candock {

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
			geometry::Point __crd;
			double __energy;
		public:
			LinkedConf(T &molecule, geometry::Point crd, double energy) 
				: __molecule(molecule), __crd(crd), __energy(energy) {}
			geometry::Point& crd() { return __crd; }
			const geometry::Point& crd() const { return __crd; }
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
		static geometry::Point::Vec greedy(const geometry::Point::Vec &initial, const double clus_rad);
		static Molib::Molecules greedy(const Molib::Molecules &initial, const Score::Score &score, Molib::Atom::Grid &gridrec, const double clus_rad);
		static Linker::Partial::Vec greedy(const Linker::Partial::Vec &initial, const Molib::Atom::Grid &gridrec, const double clus_rad);
	};
		
};

}

#endif
