#ifndef DOCKFRAGMENTS_H
#define DOCKFRAGMENTS_H

#include "programstep.hpp"

#include "findcentroids.hpp"
#include "fragmentligands.hpp"

#include "opts_candock.hpp"
#include "score/score.hpp"
#include "docker/dock.hpp"
#include "pdbreader/nrset.hpp"

namespace Program {

	class DockFragments : public ProgramStep
	{
		const FindCentroids& __found_centroids;
		const FragmentLigands& __fragmented_ligands;
		
		const Molib::Score& __score;
		const Molib::Atom::Grid& __gridrec;
		
		const std::string& __name;

		Molib::NRset __all_seeds;

		void __dock_fragment(int start, const Docker::Gpoints& gpoints, const Docker::Gpoints& gpoints0, const CmdLnOpts& cmdl);
	protected:
		virtual bool __can_read_from_files(const CmdLnOpts& cmdl);
		virtual void __read_from_files(const CmdLnOpts& cmdl);
		virtual void __continue_from_prev(const CmdLnOpts& cmdl);
		
	public:
		DockFragments ( const FindCentroids& found_centroids,
						const FragmentLigands& fragmented_ligands,
						const Molib::Score& score,
						const Molib::Atom::Grid& gridrec,
						const std::string& name
 						) :
						__found_centroids(found_centroids),
						__fragmented_ligands(fragmented_ligands),
						__score(score),
						__gridrec(gridrec),
						__name(name)
						{}

	};

}

#endif // DOCKFRAGMENTSSTEP_H
