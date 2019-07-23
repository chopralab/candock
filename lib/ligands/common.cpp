/* This is common.cpp and is part of CANDOCK
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

#include "candock/ligands/common.hpp"
#include "candock/external/json/json.h"
#include "statchem/helper/help.hpp"
#include "statchem/molib/molecule.hpp"
#include "statchem/parser/fileparser.hpp"

using namespace std;

namespace candock {

ostream& operator<<(ostream& stream, const ResidueSet& residues) {
    for (auto& presidue : residues) stream << *presidue << endl;
    return stream;
}

namespace common_ligands {

using namespace statchem;

ResMap json_to_map(const Json::Value& aligned_residues) {
    ResMap a;
    // convert Json aligned residues to a set of res_tuple2
    for (auto& rpair : aligned_residues) {
        vector<string> tok = help::ssplit(rpair["a"].asString(), "=");
        vector<string> tok1 = help::ssplit(tok[0], ":");
        vector<string> tok2 = help::ssplit(tok[1], ":");
        a.insert(make_pair(molib::Residue::res_tuple2(tok1[2].at(0), tok1[0],
                                                      stoi(tok1[1]), ' '),
                           molib::Residue::res_tuple2(
                               tok2[2].at(0), tok2[0], stoi(tok2[1]),
                               ' ')));  // {"a":"Y:212:A=Y:212:B","c":"x"}
        dbgmsg("A " << tok2[2].at(0) << ":" << tok2[0] << ":" << stoi(tok2[1]));
    }
    return a;
}
pair<ResSet, ResSet> json_to_set(Json::Value aligned_residues) {
    ResSet a, b;
    // convert Json aligned residues to a set of res_tuple2
    for (auto& rpair : aligned_residues) {
        vector<string> tok = help::ssplit(rpair["a"].asString(), "=");
        vector<string> tok1 = help::ssplit(tok[0], ":");
        vector<string> tok2 = help::ssplit(tok[1], ":");
        a.insert(molib::Residue::res_tuple2(
            tok1[2].at(0), tok1[0], stoi(tok1[1]),
            ' '));  // {"a":"Y:212:A=Y:212:B","c":"x"}
        b.insert(molib::Residue::res_tuple2(
            tok2[2].at(0), tok2[0], stoi(tok2[1]),
            ' '));  // {"a":"Y:212:A=Y:212:B","c":"x"}
        dbgmsg("A " << tok2[2].at(0) << ":" << tok2[0] << ":" << stoi(tok2[1]));
    }
    return {a, b};
}
ResMap json_to_map_reverse(Json::Value aligned_residues) {
    ResMap result;
    // convert Json aligned residues to a set of res_tuple2
    for (auto& rpair : aligned_residues) {
        vector<string> tok = help::ssplit(rpair["a"].asString(), "=");
        vector<string> tok1 = help::ssplit(tok[0], ":");
        vector<string> tok2 = help::ssplit(tok[1], ":");
        result.insert({molib::Residue::res_tuple2(tok2[2].at(0), tok2[0],
                                                  stoi(tok2[1]), ' '),
                       molib::Residue::res_tuple2(
                           tok1[2].at(0), tok1[0], stoi(tok1[1]),
                           ' ')});  // {"a":"Y:212:A=Y:212:B","c":"x"}
        dbgmsg("A " << tok2[2].at(0) << ":" << tok2[0] << ":" << stoi(tok2[1]));
    }
    //~ return {a, b};
    return result;
}
}
}
