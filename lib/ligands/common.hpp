#ifndef COMMON_HPP
#define COMMON_HPP

namespace Molib {
	class Residue;
	class Molecules;
}
typedef set<pair<Molib::Residue::res_tuple2, Molib::Residue::res_tuple2>> ResPairSet;
typedef map<Molib::Residue::res_tuple2, Molib::Residue::res_tuple2> ResMap;
typedef set<Molib::Residue::res_tuple2> ResSet;
typedef vector<Molib::Residue::res_tuple2> ResVec;
typedef set<Molib::Residue*> ResidueSet;

ostream& operator<< (ostream& stream, const ResidueSet& residues);

namespace common_ligands {
	void output_file(Molib::Molecules &mols, const string &filename);
	ResMap json_to_map(const Json::Value &aligned_residues);
	pair<ResSet, ResSet> json_to_set(Json::Value aligned_residues);
	ResMap json_to_map_reverse(Json::Value aligned_residues);
};

#endif
