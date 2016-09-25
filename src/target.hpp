#ifndef TARGET_H
#define TARGET_H

#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "docker/gpoints.hpp"
#include "pdbreader/molecules.hpp"
#include "opts_candock.hpp"

#include "findcentroids.hpp"
#include "dockfragments.hpp"

#include <string>

namespace Program {

	// TODO: Implement as a templated_map<Molecule,Target,Target> ????? Or as a Molecules?????

	class Target {

		// FIXME: There's a better to design this, but this works for *now*
		struct DockedReceptor {
			const Molib::Molecule& protein;
			std::unique_ptr<Molib::Atom::Grid> gridrec;
			std::unique_ptr<FindCentroids>     centroids;
			std::unique_ptr<DockFragments>     prepseeds;
		};

		Molib::Molecules            __receptors;
		std::vector<DockedReceptor> __preprecs;
	public:
		Target(const CmdLnOpts& cmdl, const std::string& input_name);

		std::set<int> get_idatm_types( const std::set<int>& previous = std::set<int>() ) { 
			return __receptors.get_idatm_types(previous);
		}

		void dock_fragments(const Molib::Score& score, const FragmentLigands& ligand_fragments, const CmdLnOpts& cmdl);
	};

}

#endif // TARGET_H
