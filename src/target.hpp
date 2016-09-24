#ifndef TARGET_H
#define TARGET_H

#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "docker/gpoints.hpp"
#include "pdbreader/molecules.hpp"
#include "opts_candock.hpp"

#include "findcentroids.hpp"

#include <string>

namespace Program {

	class Target {
		Molib::Molecules   __receptors;
		Molib::Atom::Grid  __gridrec;
		std::vector<FindCentroids> __centroids;
	public:
		Target(const CmdLnOpts& cmdl);
	};

}

#endif // TARGET_H
