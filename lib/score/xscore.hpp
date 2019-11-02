/* This is xscore.hpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use.
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

#ifndef XSCORE_H
#define XSCORE_H

#include <vector>
#include <array>
#include "molib/molecule.hpp"

namespace Molib {

std::array<double, 5> vina_xscore(const Molib::Atom::Grid& gridrec,
                                  const Molib::Atom::Vec& atoms);

}

#endif
