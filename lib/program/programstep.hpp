#ifndef PROGRAMSTEP_H
#define PROGRAMSTEP_H

#include "cmdlnopts.hpp"
#include <boost/filesystem.hpp>

namespace Program {

	// TODO: Possibly introduce an iterator function to iterator over results
	class ProgramStep {
	protected:
		virtual bool __can_read_from_files(const CmdLnOpts& cmdl) = 0;
		virtual void __read_from_files(const CmdLnOpts& cmdl) = 0;
		virtual void __continue_from_prev(const CmdLnOpts& cmdl) = 0;

	public:

		void run_step(const CmdLnOpts& cmdl) {
			if ( __can_read_from_files(cmdl) ) {
				__read_from_files(cmdl);
			} else {
				__continue_from_prev(cmdl);
			}
		}

	};

}

#endif
