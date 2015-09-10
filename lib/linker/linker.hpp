#ifndef LINKER_H
#define LINKER_H
#include "helper/debug.hpp"
#include "pdbreader/it.hpp"
#include "fragmenter/fragmenter.hpp"
#include "geom3d/coordinate.hpp"
#include "pdbreader/internal.hpp"
#include "pdbreader/grid.hpp"
#include <tuple>
#include <functional>
#include "segment.hpp"
#include "seed.hpp"
#include "partial.hpp"
#include "helper/array2d.hpp"

namespace Geom3D {
	class Quaternion;
};

namespace Molib {
	class Atom;
	class Molecule;
	class Molecules;
	class NRset;
	class Score;
};

namespace OMMIface {
	class Modeler;
};

namespace Linker {
	class State;
	
	class DockedConformation {
	public:
		typedef vector<DockedConformation> Vec;
	private:
		unique_ptr<Molib::Molecule> __ligand;
		unique_ptr<Molib::Molecule> __receptor;
		double __energy;
	public:
		DockedConformation() : __ligand(nullptr), __receptor(nullptr), __energy(0) {}
		DockedConformation(Molib::Molecule ligand, Molib::Molecule receptor, double energy) 
			: __ligand(unique_ptr<Molib::Molecule>(new Molib::Molecule(ligand))), 
			__receptor(unique_ptr<Molib::Molecule>(new Molib::Molecule(receptor))), __energy(energy) {}
		Molib::Molecule& get_ligand() { return *__ligand; }
		const Molib::Molecule& get_ligand() const { return *__ligand; }
		Molib::Molecule& get_receptor() { return *__receptor; }
		const Molib::Molecule& get_receptor() const { return *__receptor; }
		double get_energy() { return __energy; }
		const double get_energy() const { return __energy; }
		friend ostream& operator<<(ostream& os, const DockedConformation &conf);
	};
		
	class Linker {
		typedef map<const Molib::Atom*, Segment*> AtomToSegment;
		typedef map<const Segment*, State*> SegStateMap;
		static const int MAX_ENERGY = 999999;
		typedef multiset<Partial, Partial::comp> PriorityQueue;

	public:
		class ConnectionError : public Error {
			const set<State::ConstPair> __fs; 
		public: 
			ConnectionError(const string &msg, const set<State::ConstPair> fs) : Error(msg), __fs(fs) {}
			const set<State::ConstPair>& get_failed_state_pairs() const { return __fs; }
		};
	private:
		Molib::Internal &__ic;
		OMMIface::Modeler &__modeler;
		const Molib::Molecule &__receptor, &__ligand;
		const Molib::NRset &__top_seeds;
		Molib::Atom::Grid &__gridrec;
		const Molib::Score &__score;
		const double __dist_cutoff, __spin_degrees,
			__tol_seed_dist, __clash_coeff, __docked_clus_rad,
			__max_allow_energy;
		const int __max_possible_conf, __link_iter;
		const bool __iterative;
		const int __max_num_possibles;
		
		Segment::Graph __segment_graph;
		Seed::Graph __seed_graph;
		set<State::ConstPair> __blacklist;
		

		double __distance(const State &start, const State &goal) const;
		State::Vec __compute_neighbors(const State &curr_state, Segment &next,
		vector<unique_ptr<State>> &states);
		bool __clashes_receptor(const State&) const;
		bool __clashes_ligand(const State &current, 
			const Partial &conformation, const State &prev) const;
		Geom3D::Point::Vec __rotate(const Geom3D::Quaternion &q, const Geom3D::Point &p1, 
			const Geom3D::Point &p2, const Geom3D::Point::Vec &crds);

		void __create_states(const Segment::Graph &segment_graph, const Molib::NRset &top_seeds);
		Array2d<bool> __find_compatible_state_pairs(const Seed::Graph &seed_graph, const int sz);
		vector<vector<State::Vec>> __grow_possibles(const map<State*, State::Set> &pos);
		Partial::Vec __generate_rigid_conformations(const Seed::Graph &seed_graph);
		DockedConformation __reconstruct(const Partial &conformation);

		DockedConformation::Vec __connect(const int segment_graph_size, const Partial::Vec &possibles);
		bool __has_blacklisted(const State::Vec &conformation, const set<State::ConstPair> &blacklist);
		bool __check_distances_to_seeds(const State &curr_state, 
			const Segment &adjacent, const SegStateMap &docked_seeds);
		DockedConformation __a_star_static(const int segment_graph_size, 
			const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter);
		DockedConformation __a_star_iterative(const int segment_graph_size, 
			const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter);
		string to_pdb(const Partial &conformation);
		pair<State*, Segment*> __find_good_neighbor(const Partial &curr_conformation, 
			const SegStateMap &docked_seeds);
		State* __is_seed(const Segment &seg, const SegStateMap &docked_seeds);
	public:
		Linker(OMMIface::Modeler &modeler, const Molib::Molecule &receptor, const Molib::Molecule &ligand, 
			const Molib::NRset &top_seeds, Molib::Atom::Grid &gridrec, 
			const Molib::Score &score, Molib::Internal &ic, const double dist_cutoff, 
			const double spin_degrees, const double tol_seed_dist, 
			const int max_possible_conf, const int link_iter, const double clash_coeff, const double docked_clus_rad,
			const double max_allow_energy, bool iterative, const int max_num_possibles) : 
			__modeler(modeler), __receptor(receptor), __ligand(ligand), 
			__top_seeds(top_seeds), __gridrec(gridrec), __score(score), __ic(ic), 
			__dist_cutoff(dist_cutoff), __spin_degrees(Geom3D::radians(spin_degrees / 2)), 
			__tol_seed_dist(tol_seed_dist), __max_possible_conf(max_possible_conf),
			__link_iter(link_iter), __clash_coeff(clash_coeff), __docked_clus_rad(docked_clus_rad),
			__max_allow_energy(max_allow_energy), __iterative(iterative), __max_num_possibles(max_num_possibles) {}
		DockedConformation::Vec connect();
		Partial::Vec init_conformations();
		DockedConformation compute_conformation(const Partial &partial);
	};
}
#endif
