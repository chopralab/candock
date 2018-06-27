#ifndef LINKER_H
#define LINKER_H
#include "candock/helper/debug.hpp"
#include "candock/molib/it.hpp"
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/geometry/coordinate.hpp"
#include "candock/molib/internal.hpp"
#include "candock/molib/grid.hpp"
#include <tuple>
#include <functional>
#include "candock/helper/array2d.hpp"
#include "candock/linker/segment.hpp"
#include "candock/linker/seed.hpp"
#include "candock/linker/partial.hpp"
#include "candock/linker/dockedconformation.hpp"

namespace geometry {
        class Quaternion;
};

namespace Molib {
        class Atom;
        class Molecule;
        class Molecules;
        class NRset;
}

namespace Score {
        class Score;
}

namespace OMMIface {
        class Modeler;
}

namespace Linker {
	class State;
	
	class Linker {
	
		class GenericLinker {
        public:
            typedef map<const Segment*, State*> SegStateMap;
            
		protected:
			typedef map<const Molib::Atom*, Segment*> AtomToSegment;
			
			static const int MAX_ENERGY = 999999;
			typedef multiset<Partial, Partial::comp> PriorityQueue;
	
			class ConnectionError : public Error {
				const set<State::ConstPair> __fs; 
			public: 
				ConnectionError(const string &msg, const set<State::ConstPair> fs) : Error(msg), __fs(fs) {}
				const set<State::ConstPair>& get_failed_state_pairs() const { return __fs; }
			};
			
			Molib::Internal __ic;
			OMMIface::Modeler &__modeler;
			const Molib::Molecule &__receptor, &__ligand;
			const Molib::NRset &__top_seeds;
			const Molib::Atom::Grid &__gridrec;
			const Score::Score &__score;
			const double __dist_cutoff, __spin_degrees,
				__tol_seed_dist, __lower_tol_seed_dist, __upper_tol_seed_dist, 
				__clash_coeff, __docked_clus_rad,
				__max_allow_energy;
			const int __max_possible_conf, __link_iter;
			const int __max_num_possibles;

			const int __max_clique_size;
	
			const int __max_iterations_final;
			const int __max_iterations_pre;

			const std::string& __platform, __precision, __accelerators;
			
			Segment::Graph __segment_graph;
			Seed::Graph __seed_graph;
			set<State::ConstPair> __blacklist;
			
	
			double __distance(const State &start, const State &goal) const;
			State::Vec __compute_neighbors(const State &curr_state, Segment &next,
			vector<unique_ptr<State>> &states);
			bool __clashes_receptor(const State&) const;
			bool __clashes_ligand(const State &current, 
				const Partial &conformation, const State &prev) const;
			geometry::Point::Vec __rotate(const geometry::Quaternion &q, const geometry::Point &p1, 
				const geometry::Point::Vec &crds);
	
			Array2d<bool> __find_compatible_state_pairs(const Seed::Graph &seed_graph, const int sz);
			vector<vector<State::Vec>> __grow_possibles(const map<State*, State::Set> &pos);
	
			bool __has_blacklisted(const State::Vec &conformation, const set<State::ConstPair> &blacklist);
			bool __check_distances_to_seeds(const State &curr_state, 
				const Segment &adjacent, const SegStateMap &docked_seeds);
	
			string to_pdb(const Partial &conformation);
			pair<State*, Segment*> __find_good_neighbor(const Partial &curr_conformation, 
				const SegStateMap &docked_seeds);
			State* __is_seed(const Segment &seg, const SegStateMap &docked_seeds);

			DockedConformation::Vec __minimize_final(DockedConformation::Vec &docked_conformations);
			void __create_states(const Segment::Graph &segment_graph, const Molib::NRset &top_seeds);
			
			virtual DockedConformation __a_star(const int segment_graph_size, const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter) = 0;
			virtual Partial::Vec __generate_rigid_conformations(const Seed::Graph &seed_graph) = 0;
			virtual DockedConformation __reconstruct(const Partial &conformation) = 0;

		public:
			GenericLinker(OMMIface::Modeler &modeler, const Molib::Molecule &receptor, 
				const Molib::Molecule &ligand, const Molib::NRset &top_seeds, 
				const Molib::Atom::Grid &gridrec, const Score::Score &score, 
				const double dist_cutoff, const double spin_degrees, 
				const double tol_seed_dist, const double lower_tol_seed_dist, 
				const double upper_tol_seed_dist, const int max_possible_conf, 
				const int link_iter, const double clash_coeff, 
				const double docked_clus_rad, const double max_allow_energy, 
				const int max_num_possibles, 
				const int max_clique_size, const int max_iterations_final, const int max_iterations_pre,
				const std::string& platform, const std::string& precision, const std::string& accelerators) : 
				__ic(ligand.get_atoms()), __modeler(modeler), __receptor(receptor), 
				__ligand(ligand), __top_seeds(top_seeds), __gridrec(gridrec), 
				__score(score), __dist_cutoff(dist_cutoff),
				__spin_degrees(geometry::radians(spin_degrees / 2)), // it needs to be divided by two due to quaternion rotation (it seems to double the angles)..
				__tol_seed_dist(tol_seed_dist), __lower_tol_seed_dist(lower_tol_seed_dist), 
				__upper_tol_seed_dist(upper_tol_seed_dist), __clash_coeff(clash_coeff),
				__docked_clus_rad(docked_clus_rad), __max_allow_energy(max_allow_energy),
				__max_possible_conf(max_possible_conf), __link_iter(link_iter),
				__max_num_possibles(max_num_possibles), 
				__max_clique_size(max_clique_size), __max_iterations_final(max_iterations_final), __max_iterations_pre(max_iterations_pre),
				__platform(platform), __precision(precision), __accelerators(accelerators)
			{}

			virtual ~GenericLinker() {}

			void init_openmm();
			Partial::Vec init_conformations();
			DockedConformation::Vec compute_conformations(const Partial::Vec &partials);
			DockedConformation::Vec compute_conformations(Partial::Vec &partials);
		};
		
		class StaticLinker : public GenericLinker {
			DockedConformation __a_star(const int segment_graph_size, const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter);
			Partial::Vec __generate_rigid_conformations(const Seed::Graph &seed_graph);
			DockedConformation __reconstruct(const Partial &conformation);
		public:
			using GenericLinker::GenericLinker;
		};
		
		class IterativeLinker : public GenericLinker {
			DockedConformation __a_star(const int segment_graph_size, const Partial &start_conformation, vector<unique_ptr<State>> &states, int iter);
			Partial::Vec __generate_rigid_conformations(const Seed::Graph &seed_graph);
			DockedConformation __reconstruct(const Partial &conformation);

		public:
			using GenericLinker::GenericLinker;
		};
		

		GenericLinker *l;
		
		
	public:
		Linker(OMMIface::Modeler &modeler, const Molib::Molecule &receptor, const Molib::Molecule &ligand,
			const Molib::NRset &top_seeds, const Molib::Atom::Grid &gridrec,
			const Score::Score &score, const bool cuda, const bool iterative, const double dist_cutoff,
			const double spin_degrees, const double tol_seed_dist,
			const double lower_tol_seed_dist, const double upper_tol_seed_dist,
			const int max_possible_conf, const int link_iter,
			const double clash_coeff, const double docked_clus_rad,
			const double max_allow_energy, const int max_num_possibles,
			const int max_clique_size, const int max_iterations_final, const int max_iterations_pre,
			const std::string& platform, const std::string& precision, const std::string& accelerators);
			
		~Linker() { delete l; }
		DockedConformation::Vec link();
	};

};
#endif
