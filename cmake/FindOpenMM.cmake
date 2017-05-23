# Find openmm
#
# Find the openmm includes and library
# 
# if you nee to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  OPENMM_INCLUDE_DIR, where to find header, etc.
#  OPENMM_LIBRARY, the libraries needed to use openmm.
#  OPENMM_FOUND, If false, do not try to use openmm.

# only look in default directories
find_path(
    OPENMM_INCLUDE_DIR 
    NAMES OpenMM.h
    PATHS ${OPEMMM_ROOT}/include /usr/local/openmm/include
    DOC "openmm include dir"
)

find_library(
    OPENMM_LIBRARY
    NAMES OpenMM
    PATHS ${OPEMMM_ROOT}/lib /usr/local/openmm/lib/
    DOC "openmm library"
)

find_library(
    OPENMM_KBPLUGIN
    NAMES KBPlugin
    PATHS ${OPEMMM_ROOT}/lib /usr/local/openmm/lib/
    DOC "KBPlugin Plugin"
)

if(CUDA_FOUND)

    find_library(
        OPENMM_LIBRARY_CUDA
        NAMES OpenMMCUDA
        PATHS ${OPEMMM_ROOT}/plugins /usr/local/openmm/plugins/
        DOC "openmm library cuda"
    )

endif(CUDA_FOUND)

set(OPENMM_INCLUDE_DIR ${OPENMM_INCLUDE_DIR} CACHE STRING INTERNAL)
set(OPENMM_LIBRARY ${OPENMM_LIBRARY} ${OPENMM_LIBRARY_CUDA} ${OPENMM_KBPLUGIN} CACHE STRING INTERNAL)

IF(OPENMM_LIBRARY AND OPENMM_INCLUDE_DIR AND OPENMM_KBPLUGIN)

SET(OPENMM_FOUND 1)

ENDIF(OPENMM_LIBRARY AND OPENMM_INCLUDE_DIR AND OPENMM_KBPLUGIN)

# ==========================================
IF(NOT OPENMM_FOUND)
    MESSAGE(FATAL_ERROR "OPENMM required, please specify it's location. Try setting OPEMMM_ROOT.")
ENDIF(NOT OPENMM_FOUND)
