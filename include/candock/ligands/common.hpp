/* This is common.hpp and is part of CANDOCK
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

#ifndef COMMON_HPP
#define COMMON_HPP

#include "statchem/molib/residue.hpp"

namespace statchem {
    namespace molib {
        class Molecules;
    }
}

namespace Json {
    class Value;
}

namespace candock {

typedef std::set<
    std::pair<statchem::molib::Residue::res_tuple2, statchem::molib::Residue::res_tuple2>>
    ResPairSet;
typedef std::map<statchem::molib::Residue::res_tuple2, statchem::molib::Residue::res_tuple2> ResMap;
typedef std::set<statchem::molib::Residue::res_tuple2> ResSet;
typedef std::vector<statchem::molib::Residue::res_tuple2> ResVec;
typedef std::set<statchem::molib::Residue*> ResidueSet;

std::ostream& operator<<(std::ostream& stream, const ResidueSet& residues);

namespace common_ligands {

void output_file(statchem::molib::Molecules& mols, const std::string& filename);
ResMap json_to_map(const Json::Value& aligned_residues);
std::pair<ResSet, ResSet> json_to_set(Json::Value aligned_residues);
ResMap json_to_map_reverse(Json::Value aligned_residues);

}
}

#endif
