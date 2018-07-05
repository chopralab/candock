#ifndef COMMON_HPP
#define COMMON_HPP

#include "candock/molib/residue.hpp"

namespace candock {

namespace molib {
	class Residue;
	class Molecules;
}
typedef set<pair<molib::Residue::res_tuple2, molib::Residue::res_tuple2>> ResPairSet;
typedef map<molib::Residue::res_tuple2, molib::Residue::res_tuple2> ResMap;
typedef set<molib::Residue::res_tuple2> ResSet;
typedef vector<molib::Residue::res_tuple2> ResVec;
typedef set<molib::Residue*> ResidueSet;

ostream& operator<< (ostream& stream, const ResidueSet& residues);

namespace common_ligands {
	void output_file(molib::Molecules &mols, const string &filename);
	ResMap json_to_map(const Json::Value &aligned_residues);
	pair<ResSet, ResSet> json_to_set(Json::Value aligned_residues);
	ResMap json_to_map_reverse(Json::Value aligned_residues);
};

}

#endif
