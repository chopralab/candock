#ifndef INTERFACE_H
#define INTERFACE_H

#include "program/candockexport.hpp"

#ifdef __cplusplus
extern "C" {
#endif

const char* CANDOCK_EXPORT initialize_receptor(const char* filename);
const char* CANDOCK_EXPORT initialize_ligand(const char* filename);
      void  CANDOCK_EXPORT initialize_scoring(const char* filename);

#ifdef __cplusplus
}
#endif


#endif