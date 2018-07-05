#ifndef GENLIG_H
#define GENLIG_H
#include <iostream>
#include "candock/molib/molecule.hpp"
using namespace std;

namespace candock {

namespace genlig {
	typedef map<int, molib::Molecules> BindingSiteClusters;
	typedef map<int, double> BindingSiteScores;
	void generate_ligands(const string &receptor_file, const string &receptor_chain_id, 
		const string &json_file, const string &bio_dir, const string &lig_code,
		const string &lig_file, const string &bsite_file);
	pair<BindingSiteClusters, BindingSiteScores> generate_binding_site_prediction(const string &json_with_ligs_file, 
		const string &bio_dir, const int num_bsites);

}

// Operator overloading for typedef types can lead to issues.
// See http://blog.mezeske.com/?p=170 for details.

namespace molib {
    ostream& operator<<(ostream& os, const map<int, molib::Molecules>& bsites);
}

	ostream& operator<<(ostream& os, const map<int, double>& bscores);
}

#endif
