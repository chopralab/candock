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

	class __declspec(dllexport) Target {

		// FIXME: There's a better to design this, but this works for *now*
		// TODO:  Consider using ProgramSteps instead of named things?
		struct __declspec(dllexport) DockedReceptor {

			DockedReceptor(Molib::Molecule& rec) : protein(rec) {}

			virtual ~DockedReceptor() {
				delete score;
				delete ffield;
				delete gridrec;
				delete centroids;
				delete prepseeds;
				delete dockedlig;
			}

			DockedReceptor(const DockedReceptor& rhs) = default;
			DockedReceptor& operator= (const DockedReceptor& rhs) = delete;

			Molib::Molecule&      protein;
			Molib::Score*         score;
			OMMIface::ForceField* ffield;
			Molib::Atom::Grid*    gridrec;
			FindCentroids*        centroids;
			DockFragments*        prepseeds;
			LinkFragments*        dockedlig;
		};

		Molib::Molecules            __receptors;
		std::vector<DockedReceptor> __preprecs;
	public:
		Target(const std::string& input_name);
		Target(const Target &) = delete;
		Target& operator=(const Target&) = delete;

		std::set<int> get_idatm_types(const std::set<int>& previous = std::set<int>());

		// TODO: Instead of named function, pass in fully initiallized ProgramSteps????????
		void find_centroids();
                void rescore_docked(const FragmentLigands& ligand_fragments);
		void dock_fragments(const FragmentLigands& ligand_fragments);
		void link_fragments();
                void make_scaffolds(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add);
		void design_ligands(FragmentLigands& ligand_fragments, const std::set<std::string>& seeds_to_add);

		static void make_objective();

		// TODO: Ideally this would be done internally.....
		std::multiset<std::string> determine_overlapping_seeds(const int max_seeds, const int number_of_occurances) const;

		static std::set<std::string> determine_non_overlapping_seeds( const Target& targets, const Target& antitargets );
	};

}

#endif // TARGET_H
