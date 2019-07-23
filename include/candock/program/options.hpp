/* This is options.hpp and is part of CANDOCK
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

#ifndef OPTIONS_H
#define OPTIONS_H

#include <memory>
#include <vector>

#include "candock/candockexport.hpp"

namespace candock {

namespace help {
class CANDOCK_EXPORT Options {
   private:
    static std::unique_ptr<Options> __current_options;

   public:
    virtual ~Options();
    static const Options* get_options();
    static void set_options(Options* opts);

    virtual const std::string& get_string_option(
        const std::string& option) const = 0;
    virtual bool get_bool_option(const std::string& option) const = 0;
    virtual int get_int_option(const std::string& option) const = 0;
    virtual double get_double_option(const std::string& option) const = 0;

    virtual const std::vector<std::string>& get_string_vector(
        const std::string& option) const = 0;

    virtual std::string program_name() const = 0;
    virtual int ncpu() const = 0;

    virtual std::string configuration_file() const = 0;
};
}

#define cmdl (*candock::help::Options::get_options())
}

#endif