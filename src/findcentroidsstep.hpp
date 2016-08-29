#ifndef FINDCENTROIDSSTEP_H
#define FINDCENTROIDSSTEP_H

#include "programstep.hpp"
#include "centro/centroids.hpp"
#include "pdbreader/molecule.hpp"

namespace Program {

	class FindCentroidsStep : public ProgramStep< Centro::Centroids >
	{
	protected:
		virtual bool __can_read_from_files(const CmdLnOpts& opts);
		virtual void __read_from_files(const CmdLnOpts& cmdl);
		virtual void __continue_from_prev(const CmdLnOpts& cmdl, const ProgramStep* prev );
		
		const Molib::Molecule& __receptor;

		Centro::Centroids __result;

	public:
		FindCentroidsStep( const Molib::Molecule& receptor ) :
			__receptor( receptor ) { }

		virtual const Centro::Centroids& get_results() const {
			return __result;
		}

	};

}

#endif
