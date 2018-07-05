#include "candock/linker/dockedconformation.hpp"
#include "candock/molib/molecule.hpp"

using namespace std;

namespace candock{
namespace Linker {	
	ostream& operator<<(ostream& os, const DockedConformation &conf)	{
		os << "start link ++++++++++++++++++++++++++++++" << endl;
		os << "ENERGY = " << conf.get_energy() << endl;
		os << "LIGAND = " << conf.get_ligand() << endl;
		os << "RECEPTOR = " << conf.get_receptor() << endl;
		os << "end link --------------------------------" << endl;
		return os;
	}

	void DockedConformation::sort(DockedConformation::Vec& v) {
		std::sort(v.begin(), v.end());
	}

}
}
