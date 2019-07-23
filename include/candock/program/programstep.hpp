/* This is programstep.hpp and is part of CANDOCK
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
        if (__can_read_from_files()) {
            __read_from_files();
        } else {
            __continue_from_prev();
        }
    }
};
}
}

#endif
