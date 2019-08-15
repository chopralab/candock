/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

#include "options.hpp"

namespace help {
        std::unique_ptr<Options> Options::__current_options(nullptr);

        Options::~Options() {
                // Prevent a double free
                __current_options.release();
        }

        const help::Options * Options::get_options() {
                return __current_options.get();
        }

        void Options::set_options(help::Options* opts) {
                __current_options.reset(opts);
        }

}
