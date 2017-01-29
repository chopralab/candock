#ifndef FINDCENTROIDS_H
#define FINDCENTROIDS_H

#include "programstep.hpp"
#include "centro/centroids.hpp"
#include "pdbreader/molecule.hpp"

namespace Program {

	class FindCentroids : public ProgramStep
	{
	protected:
		virtual bool __can_read_from_files();
		virtual void __read_from_files();
		virtual void __continue_from_prev();
		
		const Molib::Molecule& __receptor;

		Centro::Centroids __result;

	public:
		FindCentroids ( const Molib::Molecule& receptor ) :
			__receptor( receptor ) { }

		const Centro::Centroids& centroids() const {
			return __result;
		}

	};

}

#endif
