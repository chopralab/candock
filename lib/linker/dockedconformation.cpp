/* This is dockedconformation.cpp and is part of CANDOCK
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

#include "candock/linker/dockedconformation.hpp"
#include "statchem/molib/molecule.hpp"

using namespace std;

namespace candock {
namespace linker {
ostream& operator<<(ostream& os, const DockedConformation& conf) {
    os << "start link ++++++++++++++++++++++++++++++" << endl;
    os << "ENERGY = " << conf.get_energy() << endl;
    os << "LIGAND = " << conf.get_ligand() << endl;
    os << "RECEPTOR = " << conf.get_receptor() << endl;
    os << "end link --------------------------------" << endl;
    return os;
}

void DockedConformation::sort(DockedConformation::Vec& v) {
    std::sort(v.begin(), v.end());
}
}
}
