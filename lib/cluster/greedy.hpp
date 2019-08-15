/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef GREEDY_CLUSTER_H
#define GREEDY_CLUSTER_H

//~ #include "molib/molecule.hpp"
#include "geom3d/geom3d.hpp"
#include "linker/partial.hpp"

namespace Molib {
	class Score;
	class Molecule;
	class Molecules;
	
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
		friend ostream& operator<<(ostream& os, const set<const LinkedConf<Molecule>*, LinkedConf<Molecule>::by_energy> &confs);
		friend ostream& operator<<(ostream& os, const set<const LinkedConf<Linker::Partial>*, LinkedConf<Linker::Partial>::by_energy> &confs);
	public:
		static Geom3D::Point::Vec greedy(const Geom3D::Point::Vec &initial, const double clus_rad);
		static Molib::Molecules greedy(const Molib::Molecules &initial, const Molib::Score &score,
			Molib::Atom::Grid &gridrec, const double clus_rad);
		static Linker::Partial::Vec greedy(const Linker::Partial::Vec &initial, const Molib::Atom::Grid &gridrec, const double clus_rad);
	};
		
};

#endif
