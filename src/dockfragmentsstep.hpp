#ifndef DOCKFRAGMENTSSTEP_H
#define DOCKFRAGMENTSSTEP_H

#include "programstep.hpp"

#include "findcentroidsstep.hpp"
#include "fragmentligandsstep.hpp"

#include "opts_candock.hpp"
#include "score/score.hpp"
#include "docker/dock.hpp"

namespace Program {

	class DockFragmentsStep : public ProgramStep
	{
		const FindCentroidsStep& __found_centroids;
		const FragmentLigandsStep& __fragmented_ligands;
		
		const Molib::Score& __score;
		const Molib::Atom::Grid& __gridrec;

		void __dock_fragment(int start, const Docker::Gpoints& gpoints, const Docker::Gpoints& gpoints0, const CmdLnOpts& cmdl);
	protected:
		virtual bool __can_read_from_files(const CmdLnOpts& cmdl);
		virtual void __read_from_files(const CmdLnOpts& cmdl);
		virtual void __continue_from_prev(const CmdLnOpts& cmdl);
		
	public:
		DockFragmentsStep( const FindCentroidsStep& found_centroids,
						   const FragmentLigandsStep& fragmented_ligands,
						   const Molib::Score& score,
						   const Molib::Atom::Grid& gridrec
 						) :
						   __found_centroids(found_centroids),
						   __fragmented_ligands(fragmented_ligands),
						   __score(score), __gridrec(gridrec)
						   {}

	};

}

#endif // DOCKFRAGMENTSSTEP_H
