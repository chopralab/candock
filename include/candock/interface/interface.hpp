#ifndef INTERFACE_H
#define INTERFACE_H

#include "candock/candockexport.hpp"
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

CANDOCK_EXPORT size_t initialize_complex(const char* filename);

CANDOCK_EXPORT size_t initialize_receptor(const char* filename);

CANDOCK_EXPORT size_t receptor_atom_count();
CANDOCK_EXPORT size_t receptor_atoms(size_t* idx, float* pos);
CANDOCK_EXPORT size_t receptor_atom_details(char* chain_ids, size_t* resi, size_t* rest, char* resn, size_t* elements);


CANDOCK_EXPORT size_t receptor_bond_count();
CANDOCK_EXPORT size_t receptor_bonds( size_t* bonds );

CANDOCK_EXPORT size_t initialize_ligand(const char* filename);

CANDOCK_EXPORT size_t ligand_atom_count();
CANDOCK_EXPORT size_t ligand_atoms(size_t* idx, float* pos);
CANDOCK_EXPORT size_t ligand_atom_details(char* chain_ids, size_t* resi, size_t* rest, size_t* elements);

CANDOCK_EXPORT size_t ligand_bond_count();
CANDOCK_EXPORT size_t ligand_bonds( size_t* bonds );
CANDOCK_EXPORT size_t ligand_get_neighbors( size_t atom_idx, size_t* neighbors );

CANDOCK_EXPORT size_t initialize_scoring(const char* obj_dir);

CANDOCK_EXPORT size_t initialize_scoring_full(const char* obj_dir,
                                              const char* ref, const char* func, const char* comp,
                                              float cutoff, float step, float scale);

CANDOCK_EXPORT size_t initialize_plugins(const char* plugin_dir);
CANDOCK_EXPORT size_t initialize_ffield(const char* data_dir, double dist_cutoff);
CANDOCK_EXPORT size_t initialize_modeler(const char* platform, const char* precision, const char* accelerators);

CANDOCK_EXPORT float calculate_score();

CANDOCK_EXPORT size_t set_positions_ligand  (const size_t* atoms, const float* positions, size_t size);
CANDOCK_EXPORT size_t set_positions_receptor(const size_t* atoms, const float* positions, size_t size);

CANDOCK_EXPORT size_t minimize_complex(size_t max_iter);

#ifdef __cplusplus
}
#endif

#endif
