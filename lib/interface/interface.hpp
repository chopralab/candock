#ifndef INTERFACE_H
#define INTERFACE_H

#include "program/candockexport.hpp"

#ifdef __cplusplus
extern "C" {
#endif

CANDOCK_EXPORT int initialize_receptor(const char* filename, char* receptor);
CANDOCK_EXPORT int initialize_ligand(const char* filename, char* ligand);
CANDOCK_EXPORT int initialize_scoring(const char* obj_dir);
CANDOCK_EXPORT int initialize_ffield(const char* data_dir);
CANDOCK_EXPORT float calculate_score();
CANDOCK_EXPORT int set_positions_ligand(const unsigned long* atoms, const float* positions, unsigned long size);
CANDOCK_EXPORT int set_positions_receptor(const unsigned long* atoms, const float* positions, unsigned long size);

#ifdef __cplusplus
}
#endif


#endif