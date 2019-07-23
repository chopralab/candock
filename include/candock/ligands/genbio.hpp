/* This is genbio.hpp and is part of CANDOCK
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

#ifndef GENBIO_HPP
#define GENBIO_HPP

#include <string>

namespace candock {

namespace genbio {
void generate_biological_assemblies(
    const std::string& models, const bool hydrogens,
    const std::string& pdb_dirname, const std::string& qpdb_file,
    const std::string& qcid, const bool neighb, const bool rnolig,
    const std::string& bio, const bool ralch, const std::string& json_file,
    const std::string& bio_file, const bool noalch,
    const std::string& mols_name, const bool rasym);
}
}

#endif
