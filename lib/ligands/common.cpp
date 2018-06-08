#include "candock/molib/molecule.hpp"
#include "candock/parser/fileparser.hpp"
#include "candock/ligands/common.hpp"
#include <json/json.h>

ostream& operator<< (ostream& stream, const ResidueSet& residues) {
	for (auto &presidue : residues)
		stream << *presidue << endl;
	return stream;
}

namespace common_ligands {
	ResMap json_to_map(const Json::Value &aligned_residues) {
		ResMap a;
		// convert Json aligned residues to a set of res_tuple2
		for (auto &rpair : aligned_residues) {
			vector<string> tok = help::ssplit(rpair["a"].asString(), "=");
			vector<string> tok1 = help::ssplit(tok[0], ":");
			vector<string> tok2 = help::ssplit(tok[1], ":");
			a.insert(make_pair(Molib::Residue::res_tuple2(tok1[2].at(0), tok1[0], stoi(tok1[1]), ' '),
								Molib::Residue::res_tuple2(tok2[2].at(0), tok2[0], stoi(tok2[1]), ' '))); // {"a":"Y:212:A=Y:212:B","c":"x"}
			dbgmsg("A " << tok2[2].at(0) << ":" << tok2[0] << ":" << stoi(tok2[1]));
		}
		return a;
	}
	pair<ResSet, ResSet> json_to_set(Json::Value aligned_residues) {
		ResSet a, b;
		// convert Json aligned residues to a set of res_tuple2
		for (auto &rpair : aligned_residues) {
			vector<string> tok = help::ssplit(rpair["a"].asString(), "=");
			vector<string> tok1 = help::ssplit(tok[0], ":");
			vector<string> tok2 = help::ssplit(tok[1], ":");
			a.insert(Molib::Residue::res_tuple2(tok1[2].at(0), tok1[0], stoi(tok1[1]), ' ')); // {"a":"Y:212:A=Y:212:B","c":"x"}
			b.insert(Molib::Residue::res_tuple2(tok2[2].at(0), tok2[0], stoi(tok2[1]), ' ')); // {"a":"Y:212:A=Y:212:B","c":"x"}
			dbgmsg("A " << tok2[2].at(0) << ":" << tok2[0] << ":" << stoi(tok2[1]));
		}
		return {a, b};
	}
	ResMap json_to_map_reverse(Json::Value aligned_residues) {
		ResMap result;
		// convert Json aligned residues to a set of res_tuple2
		for (auto &rpair : aligned_residues) {
			vector<string> tok = help::ssplit(rpair["a"].asString(), "=");
			vector<string> tok1 = help::ssplit(tok[0], ":");
			vector<string> tok2 = help::ssplit(tok[1], ":");
			result.insert({Molib::Residue::res_tuple2(tok2[2].at(0), tok2[0], stoi(tok2[1]), ' '), 
				Molib::Residue::res_tuple2(tok1[2].at(0), tok1[0], stoi(tok1[1]), ' ')}); // {"a":"Y:212:A=Y:212:B","c":"x"}
			dbgmsg("A " << tok2[2].at(0) << ":" << tok2[0] << ":" << stoi(tok2[1]));
		}
		//~ return {a, b};
		return result;
	}
}
