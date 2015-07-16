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
#include "helper/array2d.hpp"

namespace Geom3D {
	class Quaternion;
}
namespace Molib {
	class Atom;
	class Residue;
	class Model;
	class Molecule;
	class Molecules;
	class Segment;
	class State;
	class Score;
	class Linker {
		typedef map<const Atom*, Segment*> AtomToSegment;
		typedef map<const Segment*, State*> SegStateMap;
		typedef Grid<Atom> MolGrid;
		static const int MAX_ENERGY = 999999;
		typedef pair<State::Vec, double> LinkEnergy;
		typedef pair<State::NVec, double> Conf;
		typedef vector<Conf> Conformations;
		struct pq_comp { bool operator()(const LinkEnergy &i, const LinkEnergy &j) 
			{ return i.second < j.second; } };
		typedef multiset<LinkEnergy, pq_comp> PriorityQueue;

		class ConnectionError : public Error {
			const set<State::ConstPair> __fs; 
		public: 
			ConnectionError(const string &msg, const set<State::ConstPair> fs) : Error(msg), __fs(fs) {}
			const set<State::ConstPair>& get_failed_state_pairs() const { return __fs; }
		};
		Internal &__ic;
		const Molecule &__ligand;
		const NRset &__top_seeds;
		MolGrid &__gridrec;
		const Score &__score;
		const double __dist_cutoff, __spin_degrees, __tol_dist,
			__tol_max_coeff, __tol_min_coeff;
		const int __max_possible_conf, __link_iter;
		State::Vec __states;


		double __distance(const State &start, const State &goal) const;
		State::Vec __compute_neighbors(const State &curr_state, Segment &next,
		vector<unique_ptr<State>> &states);
		bool __clashes_receptor(const State&) const;
		bool __clashes_ligand(const State &current, 
			const LinkEnergy &conformation, const State &prev) const;
		AtomToCrd __rotate(const Geom3D::Quaternion&, const Geom3D::Point&, const Geom3D::Point&, const AtomToCrd&);

		void __create_states(const Segment::Graph &segment_graph, const NRset &top_seeds);
		Array2d<bool> __find_compatible_state_pairs(const Seed::Graph &seed_graph);
		vector<vector<State::Vec>> __grow_possibles(const map<State*, State::Set> &pos);
		vector<LinkEnergy> __generate_rigid_conformations(const Seed::Graph &seed_graph);
		Segment::Paths __find_paths(const Segment::Graph &segment_graph);
		void __set_branching_rules(const Segment::Paths &paths);
		Molecules __reconstruct(const Conformations &conformations);
		void __init_max_linker_length(const Segment::Paths &paths);
		void __compute_max_linker_length(Segment::Graph::Path &path);

		Conformations __connect(const int segment_graph_size, const vector<LinkEnergy> &possibles);
		bool __has_blacklisted(const State::Vec &conformation, const set<State::ConstPair> &blacklist);
		bool __link_adjacent(const Segment::Graph::Path &path);
		bool __check_distances_to_seeds(const State &curr_state, 
			const Segment &adjacent, const SegStateMap &docked_seeds);
		Conf __a_star(const int segment_graph_size, 
			const LinkEnergy &start_conformation, int iter);
		string to_pdb(const LinkEnergy &conformation);
		pair<State*, Segment*> __find_good_neighbor(const LinkEnergy &curr_conformation, 
			const SegStateMap &docked_seeds);
		State* __is_seed(const Segment &seg, const SegStateMap &docked_seeds);
	public:
		Linker(const Molecule &ligand, const NRset &top_seeds, MolGrid &gridrec, 
			const Score &score, Internal &ic, const double dist_cutoff, 
			const double spin_degrees, const double tol_dist, const double tol_max_coeff,
			const double tol_min_coeff, const int max_possible_conf,
			const int link_iter) : __ligand(ligand), 
			__top_seeds(top_seeds), __gridrec(gridrec), __score(score), __ic(ic), 
			__dist_cutoff(dist_cutoff), 
			__spin_degrees(Geom3D::radians(spin_degrees / 2)), 
			__tol_dist(tol_dist), __tol_max_coeff(tol_max_coeff), 
			__tol_min_coeff(tol_min_coeff), __max_possible_conf(max_possible_conf),
			__link_iter(link_iter) {}
		Molecules connect();
		friend ostream& operator<<(ostream& os, const Conf &conf);
		friend ostream& operator<<(ostream& os, const LinkEnergy &le);
		friend ostream& operator<<(ostream& os, const vector<LinkEnergy> &vec_le);
	};
}
#endif
