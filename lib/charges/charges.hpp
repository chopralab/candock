#ifndef CHARGES_H
#define CHARGES_H
#include <string>
#include <map>
#include <set>
#include <vector>
#include "geom3d/coordinate.hpp"
#include "pdbreader/molecule.hpp"
#include <stdlib.h>
using namespace std;
namespace Molib {
	class Charges {
		void __apply_am1_bcc(Molecule &molecule) {
			/*
			 * from charge.c
			 */
			wac("ANTECHAMBER_AM1BCC_PRE.AC", atomnum, atom, bondnum, bond, *cinfo,
				*minfo);
			copied_size = build_exe_path(tmpchar1, "am1bcc", sizeof tmpchar1, 1);
			strncat(tmpchar1, " -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC"
				" -f ac -p " , sizeof tmpchar1 - copied_size );
				build_dat_path(tmpchar2, "BCCPARM.DAT", sizeof tmpchar2, 1);
		        strncat(tmpchar1, tmpchar2, sizeof(tmpchar1) - strlen(tmpchar1) -1);
			strcat(tmpchar1, " -s ");
			sprintf(tmpchar2, "%d", (*cinfo).intstatus);
			strcat(tmpchar1, tmpchar2);
			if ((*cinfo).prediction_index == -1)
				(*cinfo).prediction_index = 0;
			sprintf(tmpchar, "%d", (*cinfo).prediction_index);
			strcat(tmpchar1, " -j ");
			strcat(tmpchar1, tmpchar);
			status = system(tmpchar1);
		}
		void __calculate_am1_charges(Molecule &molecule) {
			/* 
			 * from charge.c
			 */
			// determine charge of molecule
			// prepare input for sqm
			// calculate Mulliken charges with sqm 
			wsqmcrt("sqm.in", atomnum, atom, *minfo);
			build_exe_path(tmpchar, "sqm", sizeof tmpchar, 1);
			strcat(tmpchar, " -O -i sqm.in -o sqm.out");
			status = system(tmpchar);
			
			rsqmcharge("sqm.out", atomnum, atom, minfo);
			wpdb_optimized("sqm.out",atomnum,atom,2); 
		}
	public:
		
	};
}
