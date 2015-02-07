#ifndef GENLIG_H
#define GENLIG_H
#include <iostream>
#include "pdbreader/molecule.hpp"
using namespace std;

namespace genlig {
	typedef map<int, Molib::Molecules> BindingSiteClusters;
	BindingSiteClusters generate_binding_site_prediction(const string &json_with_ligs_file, 
		const string &bio_dir, const int num_bsites);
}

ostream& operator<<(ostream& os, const genlig::BindingSiteClusters& bsites);

#endif
