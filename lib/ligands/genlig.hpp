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

#ifndef GENLIG_H
#define GENLIG_H
#include <iostream>
#include "molib/molecule.hpp"
using namespace std;

namespace genlig {
	typedef map<int, Molib::Molecules> BindingSiteClusters;
	typedef map<int, double> BindingSiteScores;
	void generate_ligands(const string &receptor_file, const string &receptor_chain_id, 
		const string &json_file, const string &bio_dir, const string &lig_code,
		const string &lig_file, const string &bsite_file);
	pair<BindingSiteClusters, BindingSiteScores> generate_binding_site_prediction(const string &json_with_ligs_file, 
		const string &bio_dir, const int num_bsites);

}

// Operator overloading for typedef types can lead to issues.
// See http://blog.mezeske.com/?p=170 for details.

namespace Molib {
    ostream& operator<<(ostream& os, const map<int, Molib::Molecules>& bsites);
}

namespace std {
	ostream& operator<<(ostream& os, const map<int, double>& bscores);
}

#endif
