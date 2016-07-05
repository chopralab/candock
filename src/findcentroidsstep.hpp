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
		virtual Centro::Centroids* __read_from_files(const CmdLnOpts& cmdl);
		virtual Centro::Centroids* __continue_from_prev(const CmdLnOpts& cmdl, const ProgramStep* prev );
		
		const Molib::Molecule& __receptor;

	public:
		FindCentroidsStep( const Molib::Molecule& receptor ) :
			__receptor( receptor ) { }

	};

}

#endif
