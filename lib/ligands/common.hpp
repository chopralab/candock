/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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
