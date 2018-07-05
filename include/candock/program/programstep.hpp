#ifndef PROGRAMSTEP_H
#define PROGRAMSTEP_H

#include "candock/program/cmdlnopts.hpp"

namespace candock {

namespace Program {

	// TODO: Possibly introduce an iterator function to iterator over results
	class CANDOCK_EXPORT ProgramStep {
	protected:
		virtual bool __can_read_from_files() = 0;
		virtual void __read_from_files() = 0;
		virtual void __continue_from_prev() = 0;

	public:

		void run_step() {
			if ( __can_read_from_files() ) {
				__read_from_files();
			} else {
				__continue_from_prev();
			}
		}

	};

}

}

#endif
