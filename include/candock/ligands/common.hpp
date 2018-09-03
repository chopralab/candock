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
