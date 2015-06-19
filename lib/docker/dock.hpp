#ifndef DOCK_H
#define DOCK_H

#include "gpoints.hpp"
#include "conformations.hpp"

namespace Molib {
	class Molecule;
	//~ class Score;
};

namespace Docker {
	
	typedef vector<pair<Gpoint*, OneConformation*>> AcceptedConformations;
	
	class Dock {
		Gpoints &__gpoints;
		Conformations &__conformations;
		Molib::Molecule &__seed;
		//~ Molib::Score &__score;
	public:
		Dock(Gpoints &gpoints, Conformations &conformations,
			//~ Molib::Molecule &seed, Molib::Score &score) : __gpoints(gpoints),
			Molib::Molecule &seed) : __gpoints(gpoints),
			//~ __conformations(conformations), __seed(seed), __score(score) {}
			__conformations(conformations), __seed(seed) {}
		//~ Molib::Molecules run();
		AcceptedConformations run();
		Molib::Molecules cluster(const AcceptedConformations &accepted);
	};

};

#endif
