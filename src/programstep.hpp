#ifndef PROGRAMSTEP_H
#define PROGRAMSTEP_H

#include "opts_candock.hpp"
#include <boost/filesystem.hpp>

namespace Program {

	template <class P>
	class ProgramStep {
	protected:
		virtual bool __can_read_from_files(const CmdLnOpts& cmdl) = 0;
		virtual void __read_from_files(const CmdLnOpts& cmdl) = 0;
		virtual void __continue_from_prev(const CmdLnOpts& cmdl, const ProgramStep* prev ) = 0;
		
		int __get_file_size( boost::filesystem::path p ) {
			if ( boost::filesystem::exists(p) &&
			     boost::filesystem::is_regular_file(p) )
			{
				return boost::filesystem::file_size(p);
			} else {
				return 0;
			}
		}
		
	public:

		void run_step(const CmdLnOpts& cmdl, const ProgramStep* prev) {
			if ( __can_read_from_files(cmdl) ) {
				__read_from_files(cmdl);
			} else {
				__continue_from_prev(cmdl, prev);
			}
		}

		virtual const P& get_results() const = 0;

	};

}

#endif
