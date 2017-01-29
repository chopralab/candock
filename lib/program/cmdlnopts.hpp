#ifndef CMDLNOPTS_H
#define CMDLNOPTS_H

#include <boost/asio/ip/host_name.hpp>
#include <boost/program_options.hpp>
#include "helper/error.hpp"
#include "helper/options.hpp"
#include "version.h"
#include <thread>
#include <vector>
#include <iostream>
#include <iomanip>

namespace Program {

	class CmdLnOpts : public help::Options {
		boost::program_options::variables_map __vm;
                
                bool __quiet;
                std::string __version;
                std::string __git_version;
                std::string __program_name;
                int __ncpu;

	public:
		CmdLnOpts() : __quiet (false), __version (""), __git_version ("") {
			std::stringstream ss;

			ss << CANDOCK_MAJOR_VERSION << "." << CANDOCK_MINOR_VERSION << "."
			   << CANDOCK_TWEAK_VERSION;

			__version = ss.str();

			std::stringstream ss2;

			ss2 << CANDOCK_GIT_REFERENCE << " on branch " << CANDOCK_GIT_MYCBRANCH;

			__git_version = ss2.str();
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

		void init (int argc, char *argv[], int opts_to_parse = ALL_OPTIONS);

		// interface

		

		const std::string& get_string_option (const std::string& option) const;
		bool        get_bool_option   (const std::string& option) const;
		int         get_int_option    (const std::string& option) const;
		double      get_double_option (const std::string& option) const;

		const std::vector<std::string>& get_string_vector (const std::string& option) const;

		bool quiet() const {
			return __quiet;
		}

                int ncpu() const {
                        return __ncpu;
                }

		void print_version() const {

			std::cout << 
			"-----------------------------------------------------------------------|\n"
			"|         CCCCC    AAA   NN   NN DDDDD    OOOOO   CCCCC  KK  KK        |\n"
			"|        CC       AAAAA  NNN  NN DD  DD  OO   OO CC      KK KK         |\n"
			"|        CC      AA   AA NN N NN DD   DD OO   OO CC      KKKK          |\n"
			"|        CC      AAAAAAA NN  NNN DD   DD OO   OO CC      KK KK         |\n"
			"|         CCCCC  AA   AA NN   NN DDDDDD   OOOO0   CCCCC  KK  KK        |\n"
			"|                                                                      |\n"
			"|                             It can-dock!                             |\n"
			"|                                                                      |\n" 
			"|              Copyright 2016  Chopra Lab (chopralab.org)              |\n"
			"|                          Purdue  University                          |\n"
			"|                                                                      |\n"
			"|                      Member of the CANDIY Suite                      |\n"
			"|______________________________________________________________________|\n"
			<< std::endl;

			std::cout << "running " << __program_name << " on hostname " << boost::asio::ip::host_name() << std::endl;
			std::cout << "version " << __version      << std::endl;
			std::cout << "buildid " << __git_version  << std::endl;

			std::cout << std::endl << "Detected support for "
			          << std::thread::hardware_concurrency() << " concurrent threads."
			          << " Using " << __ncpu << " threads." << endl;

		}

		friend std::ostream &operator<< (std::ostream &stream, const CmdLnOpts &cmdl_) {
			
			cmdl_.print_version();

			for ( const auto& a : cmdl_.__vm ) {
				stream << std::setw(22)<< a.first << " = ";
				if        ( auto v = boost::any_cast<std::string>(&a.second.value()) ) {
					stream << std::setw(47) << *v;
				} else if ( auto v = boost::any_cast<int>(&a.second.value()) ) {
					stream << std::setw(47) << *v;
				} else if ( auto v = boost::any_cast<double>(&a.second.value()) ) {
					stream << std::setw(47) << *v;
				} else if ( auto v = boost::any_cast<bool>(&a.second.value()) ) {
					stream << std::setw(47) << *v;
				} else if ( auto v = boost::any_cast<std::vector<std::string>>(&a.second.value()) ) {
					std::string combination("");
					for ( const auto& s : *v ) {
						combination += s;
						combination += ", ";
					}
					stream << std::setw(47) << combination;
				}
				
				stream << std::endl;
			}

			return stream;
		}

	};
}

#endif
