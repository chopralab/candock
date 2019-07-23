/* This is cmdlnopts.hpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef CMDLNOPTS_H
#define CMDLNOPTS_H

#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>
#include "statchem/helper/error.hpp"
#include "candock/program/options.hpp"

#include "candock/candockexport.hpp"

namespace candock {

namespace Program {

class CANDOCK_EXPORT CmdLnOpts : public help::Options {
    boost::program_options::variables_map __vm;

    bool __quiet;
    std::string __program_name;
    int __ncpu;

    void __init(int argc, char* argv[], int opts_to_parse = ALL_OPTIONS);

   public:
    CmdLnOpts(int argc, char* argv[], int opts_to_parse = ALL_OPTIONS)
        : __quiet(false) {
        __init(argc, argv, opts_to_parse);
    }

    enum CMDLN_OPTS_GROUPS {
        STARTING = 1 << 0,
        PROBIS = 1 << 1,
        LIG_FRAMGENT = 1 << 2,
        FRAG_DOCKING = 1 << 3,
        SCORING = 1 << 4,
        FORCE_FIELD = 1 << 5,
        LINKING = 1 << 6,
        DESIGN = 1 << 7,
        ALL_OPTIONS = 0xFFFF
    };

    const std::string& get_string_option(const std::string& option) const;
    bool get_bool_option(const std::string& option) const;
    int get_int_option(const std::string& option) const;
    double get_double_option(const std::string& option) const;

    const std::vector<std::string>& get_string_vector(
        const std::string& option) const;

    bool quiet() const { return __quiet; }

    std::string program_name() const { return __program_name; }

    int ncpu() const { return __ncpu; }

    std::string configuration_file() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    friend CANDOCK_EXPORT std::ostream& operator<<(std::ostream& stream,
                                                   const CmdLnOpts& cmdl_);
};
}
}

#endif
