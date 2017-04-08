#ifndef INTERFACE_H
#define INTERFACE_H

#include "program/candockexport.hpp"

#ifdef __cplusplus
extern "C" {
#endif

CANDOCK_EXPORT const char* initialize_receptor(const char* filename);
CANDOCK_EXPORT const char* initialize_ligand(const char* filename);
CANDOCK_EXPORT       void  initialize_scoring(const char* obj_dir);
CANDOCK_EXPORT       void  initialize_ffield(const char* data_dir);
CANDOCK_EXPORT       float calculate_score();
CANDOCK_EXPORT       void  set_positions_ligand(const unsigned long* atoms, const float* positions, unsigned long size);

#ifdef __cplusplus
}
#endif


#endif