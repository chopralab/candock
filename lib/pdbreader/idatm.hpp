#ifndef IDATM_H
#define IDATM_H

namespace Molib {
	class Molecule;
	class Molecules;
};

namespace Idatm {

	void compute_idatm_type(Molib::Molecule &molecule);
	void compute_idatm_type(Molib::Molecules &mols);
};
#endif
