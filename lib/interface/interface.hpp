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

#ifndef INTERFACE_H
#define INTERFACE_H

#include "program/candockexport.hpp"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

CANDOCK_EXPORT const char* cd_get_error();

typedef enum {
        SINGLE_BOND = 1 << 0,
        DOUBLE_BOND = 1 << 1,
        TRIPLE_BOND = 1 << 2,
        INRING_BOND = 1 << 3,
        ROTATE_BOND = 1 << 4,
        AROMAT_BOND = 1 << 5,
} cd_bond_type;

/*
 * All functions return 0 upon failure and a non-zero number upon success
 */

CANDOCK_EXPORT size_t initialize_receptor(const char* filename);
CANDOCK_EXPORT size_t receptor_atom_count();
CANDOCK_EXPORT size_t receptor_atoms(size_t* idx, float* pos);
CANDOCK_EXPORT size_t receptor_bond_count();
CANDOCK_EXPORT size_t receptor_bonds( size_t* bonds );

CANDOCK_EXPORT size_t initialize_ligand(const char* filename);
CANDOCK_EXPORT size_t ligand_atom_count();
CANDOCK_EXPORT size_t ligand_atoms(size_t* idx, float* pos);
CANDOCK_EXPORT size_t ligand_bond_count();
CANDOCK_EXPORT size_t ligand_bonds( size_t* bonds );

CANDOCK_EXPORT size_t initialize_scoring(const char* obj_dir);
CANDOCK_EXPORT size_t initialize_plugins(const char* plugin_dir);
CANDOCK_EXPORT size_t initialize_ffield(const char* data_dir);

CANDOCK_EXPORT float calculate_score();

CANDOCK_EXPORT size_t set_positions_ligand  (const size_t* atoms, const float* positions, size_t size);
CANDOCK_EXPORT size_t set_positions_receptor(const size_t* atoms, const float* positions, size_t size);

CANDOCK_EXPORT size_t minimize_complex(size_t max_iter, size_t update_freq);

#ifdef __cplusplus
}
#endif

#endif
