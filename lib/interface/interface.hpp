#ifndef INTERFACE_H
#define INTERFACE_H

#include "program/candockexport.hpp"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * All functions return 0 upon failure and a non-zero number upon success
 */

CANDOCK_EXPORT size_t initialize_receptor(const char* filename);
CANDOCK_EXPORT size_t receptor_atom_count();
CANDOCK_EXPORT size_t receptor_atoms(size_t* idx, float* pos);

CANDOCK_EXPORT size_t receptor_string_size();
CANDOCK_EXPORT size_t copy_receptor_string(char* buffer);

CANDOCK_EXPORT size_t initialize_ligand(const char* filename);
CANDOCK_EXPORT size_t ligand_atom_count();
CANDOCK_EXPORT size_t ligand_atoms(size_t* idx, float* pos);

CANDOCK_EXPORT size_t ligand_string_size();
CANDOCK_EXPORT size_t copy_ligand_string(char* buffer);

CANDOCK_EXPORT size_t initialize_scoring(const char* obj_dir);
CANDOCK_EXPORT size_t initialize_ffield(const char* data_dir);

CANDOCK_EXPORT float calculate_score();

CANDOCK_EXPORT size_t set_positions_ligand  (const size_t* atoms, const float* positions, size_t size);
CANDOCK_EXPORT size_t set_positions_receptor(const size_t* atoms, const float* positions, size_t size);

CANDOCK_EXPORT size_t minimize_complex(size_t max_iter, size_t update_freq);

#ifdef __cplusplus
}
#endif

#endif
