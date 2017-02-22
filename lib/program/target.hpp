#ifndef TARGET_H
#define TARGET_H

#include "findcentroids.hpp"
#include "dockfragments.hpp"
#include "linkfragments.hpp"
#include "design/design.hpp"

#include "docker/gpoints.hpp"
#include "pdbreader/molecules.hpp"

#include <string>

namespace Program {

	// TODO: Implement as a templated_map<Molecule,Target,Target> ????? Or as a Molecules?????

	class Target {

		// FIXME: There's a better to design this, but this works for *now*
		// TODO:  Consider using ProgramSteps instead of named things?
		struct DockedReceptor {
                        DockedReceptor( Molib::Molecule& my_protein ) : 
                                protein(my_protein),
                                score(nullptr), ffield(nullptr), gridrec(nullptr),
                                centroids(nullptr), prepseeds(nullptr), dockedlig(nullptr)
                                {}
			Molib::Molecule& protein;
			std::unique_ptr<Molib::Score>         score;
			std::unique_ptr<OMMIface::ForceField> ffield;
			std::unique_ptr<Molib::Atom::Grid>    gridrec;
			std::unique_ptr<FindCentroids>        centroids;
			std::unique_ptr<DockFragments>        prepseeds;
			std::unique_ptr<LinkFragments>        dockedlig;
		};

		Molib::Molecules            __receptors;
		std::vector<DockedReceptor> __preprecs;
	public:
		Target(const std::string& input_name);

		std::set<int> get_idatm_types( const std::set<int>& previous = std::set<int>() ) { 
			return __receptors.get_idatm_types(previous);
		}

		// TODO: Instead of named function, pass in fully initiallized ProgramSteps????????
		void find_centroids();
		void dock_fragments(const FragmentLigands& ligand_fragments);
		void link_fragments();
		void design_ligands(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add);

		// TODO: Ideally this would be done internally.....
		std::multiset<std::string> determine_overlapping_seeds(const int max_seeds, const int number_of_occurances) const;

		static std::set<std::string> determine_non_overlapping_seeds( const Target& targets, const Target& antitargets );
	};

}

#endif // TARGET_H
