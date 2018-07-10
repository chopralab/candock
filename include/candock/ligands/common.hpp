#ifndef COMMON_HPP
#define COMMON_HPP

#include "candock/molib/residue.hpp"

namespace candock {

namespace molib {
	class Residue;
	class Molecules;
}
typedef std::set<std::pair<molib::Residue::res_tuple2, molib::Residue::res_tuple2>> ResPairSet;
typedef std::map<molib::Residue::res_tuple2, molib::Residue::res_tuple2> ResMap;
typedef std::set<molib::Residue::res_tuple2> ResSet;
typedef std::vector<molib::Residue::res_tuple2> ResVec;
typedef std::set<molib::Residue*> ResidueSet;

std::ostream& operator<< (std::ostream& stream, const ResidueSet& residues);

namespace common_ligands {
	void output_file(molib::Molecules &mols, const std::string &filename);
	ResMap json_to_map(const Json::Value &aligned_residues);
	std::pair<ResSet, ResSet> json_to_set(Json::Value aligned_residues);
	ResMap json_to_map_reverse(Json::Value aligned_residues);
};

}

#endif
