#ifndef INTERFACE_H
#define INTERFACE_H

#include "program/candockexport.hpp"

#ifdef __cplusplus
extern "C" {
#endif

const char* CANDOCK_EXPORT initialize_receptor(const char* filename);
const char* CANDOCK_EXPORT initialize_ligand(const char* filename);
      void  CANDOCK_EXPORT initialize_scoring(const char* obj_dir);
      void  CANDOCK_EXPORT initialize_ffield(const char* data_dir);
      float CANDOCK_EXPORT calculate_score();
      void  CANDOCK_EXPORT set_positions_ligand(const unsigned long* atoms, const float* positions, unsigned long size);

#ifdef __cplusplus
}
#endif


#endif