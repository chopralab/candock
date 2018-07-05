#ifndef GREEDY_CLUSTER_H
#define GREEDY_CLUSTER_H

#include "candock/geometry/geometry.hpp"
#include "candock/linker/partial.hpp"
#include "candock/score/score.hpp"

namespace candock {

namespace score {
        class Score;
}

namespace molib {
        class Molecule;
        class Molecules;
}

namespace cluster {
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
		friend ostream& operator<<(ostream& os, const set<const LinkedConf<molib::Molecule>*, LinkedConf<molib::Molecule>::by_energy> &confs);
		friend ostream& operator<<(ostream& os, const set<const LinkedConf<linker::Partial>*, LinkedConf<linker::Partial>::by_energy> &confs);
	public:
		static geometry::Point::Vec greedy(const geometry::Point::Vec &initial, const double clus_rad);
		static molib::Molecules greedy(const molib::Molecules &initial, const score::Score &score, molib::Atom::Grid &gridrec, const double clus_rad);
		static linker::Partial::Vec greedy(const linker::Partial::Vec &initial, const molib::Atom::Grid &gridrec, const double clus_rad);
	};
		
};

}

#endif
