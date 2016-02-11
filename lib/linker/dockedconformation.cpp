#include "dockedconformation.hpp"
#include "pdbreader/molecule.hpp"

using namespace std;

namespace Linker {	
	ostream& operator<<(ostream& os, const DockedConformation &conf)	{
		os << "start link ++++++++++++++++++++++++++++++" << endl;
		os << "ENERGY = " << conf.get_energy() << endl;
		os << "LIGAND = " << conf.get_ligand() << endl;
		os << "RECEPTOR = " << conf.get_receptor() << endl;
		os << "end link --------------------------------" << endl;
		return os;
	}
};
