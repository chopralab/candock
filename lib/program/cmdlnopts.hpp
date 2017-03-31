#ifndef CMDLNOPTS_H
#define CMDLNOPTS_H

#include <boost/program_options.hpp>
#include "helper/error.hpp"
#include "options.hpp"
#include <thread>
#include <vector>
#include <iostream>
#include <iomanip>

#include "candockexport.hpp"

namespace Program {

        class CANDOCK_EXPORT CmdLnOpts : public help::Options {
                boost::program_options::variables_map __vm;
                
                bool __quiet;
                std::string __program_name;
                int __ncpu;
                
                void __init (int argc, char *argv[], int opts_to_parse = ALL_OPTIONS);

        public:
                CmdLnOpts(int argc, char *argv[], int opts_to_parse = ALL_OPTIONS) : __quiet (false) {
                        __init(argc, argv, opts_to_parse);
                }

		enum CMDLN_OPTS_GROUPS {
			STARTING     = 1 << 0,
			PROBIS       = 1 << 1,
			LIG_FRAMGENT = 1 << 2,
			FRAG_DOCKING = 1 << 3,
			SCORING      = 1 << 4,
			FORCE_FIELD  = 1 << 5,
			LINKING      = 1 << 6,
			DESIGN       = 1 << 7,
			ALL_OPTIONS  = 0xFFFF
		};

		const std::string& get_string_option (const std::string& option) const;
		bool        get_bool_option   (const std::string& option) const;
		int         get_int_option    (const std::string& option) const;
		double      get_double_option (const std::string& option) const;

		const std::vector<std::string>& get_string_vector (const std::string& option) const;

                bool quiet() const {
                        return __quiet;
                }

                std::string program_name() const {
                        return __program_name;
                }

                int ncpu() const {
                        return __ncpu;
                }
                
                std::string configuration_file() const {
                        stringstream ss;
                        ss << *this;
                        return ss.str();
                }

                friend CANDOCK_EXPORT std::ostream &operator<< (std::ostream &stream, const CmdLnOpts &cmdl_);

	};
}

#endif
