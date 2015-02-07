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
		typedef Glib::Graph<Segment> SegGraph;
		typedef Glib::Graph<Seed> SeedGraph;
		typedef map<const Atom*, Segment*> AtomToSegment;
		typedef map<ConstSegPair, SegGraph::Path> Paths;
		typedef map<const Segment*, State*> SegStateMap;
		typedef Grid<Atom> MolGrid;
		static const int MAX_ENERGY = 999999;
		typedef pair<StateVec, double> LinkEnergy;
		typedef pair<NStateVec, double> Conf;
		typedef vector<Conf> Conformations;
		struct pq_comp { bool operator()(const LinkEnergy &i, const LinkEnergy &j) 
			{ return i.second < j.second; } };
		typedef multiset<LinkEnergy, pq_comp> PriorityQueue;

		class ConnectionError : public Error {
			const set<ConstStatePair> __fs; 
		public: 
			ConnectionError(const string &msg, const set<ConstStatePair> fs) : Error(msg), __fs(fs) {}
			const set<ConstStatePair>& get_failed_state_pairs() const { return __fs; }
		};
		Internal &__ic;
		const Molecule &__ligand;
		const NRset &__top_seeds;
		MolGrid &__gridrec;
		const Score &__score;
		const double __dist_cutoff, __spin_degrees, __tol_dist,
			__tol_max_coeff, __tol_min_coeff;
		const int __max_possible_conf, __link_iter;
		double __distance(const State &start, const State &goal) const;
		StateVec __compute_neighbors(const State &curr_state, Segment &next,
		vector<unique_ptr<State>> &states);
		bool __clashes_receptor(const State&) const;
		bool __clashes_ligand(const State &current, 
			const LinkEnergy &conformation, const State &prev) const;
		AtomToCrd __rotate(const Geom3D::Quaternion&, const Geom3D::Point&, const Geom3D::Point&, const AtomToCrd&);

		void __create_states(const SegGraph &segment_graph, const NRset &top_seeds);
		map<State*, StateSet> __find_compatible_state_pairs(const SeedGraph &seed_graph);
		vector<vector<StateVec>> __grow_possibles(const map<State*, StateSet> &pos);
		vector<LinkEnergy> __find_possible_states(const SeedGraph &seed_graph);
		Paths __find_paths(const SegGraph &segment_graph);
		void __set_branching_rules(const Paths &paths);
		Molecules __reconstruct(const Conformations &conformations);
		void __init_max_linker_length(const Paths &paths);
		void __compute_max_linker_length(SegGraph::Path &path);

		SegGraph __create_segment_graph(const Molecule &molecule);
		SeedGraph __create_seed_graph(const SegGraph &segment_graph, const Paths &paths);

		Conformations __connect(const int segment_graph_size, const vector<LinkEnergy> &possibles);
		bool __has_blacklisted(const StateVec &conformation, const set<ConstStatePair> &blacklist);
		bool __link_adjacent(const SegGraph::Path &path);
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
	};
}
#endif
