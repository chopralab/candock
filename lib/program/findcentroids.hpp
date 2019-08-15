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

#ifndef FINDCENTROIDS_H
#define FINDCENTROIDS_H

#include "programstep.hpp"
#include "centro/centroids.hpp"
#include "molib/molecule.hpp"

namespace Program {

	class FindCentroids : public ProgramStep
	{
	protected:
		virtual bool __can_read_from_files();
		virtual void __read_from_files();
		virtual void __continue_from_prev();
		
		const Molib::Molecule& __receptor;
                const std::string&     __filename;

		Centro::Centroids __result;
                std::string __centroid_file;

	public:
		FindCentroids ( const Molib::Molecule& receptor, const std::string& filename );
                virtual ~FindCentroids() {}

		const Centro::Centroids& centroids() const {
			return __result;
		}

	};

}

#endif
