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
    PATHS /usr/local/openmm/include
    DOC "openmm include dir"
)


if(CYGWIN)
    find_library(
        OPENMM_LIBRARY
        NAMES libOpenMM.dll.a
        PATHS /usr/local/openmm/lib
        DOC "openmm library"
    )
    find_library(
        OPENMM_KBPLUGIN
        NAMES libKBPlugin.dll.a
        PATHS /usr/local/openmm/lib/
        DOC "KBPlugin Plugin"
    )
else(CYGWIN)
    find_library(
        OPENMM_LIBRARY
        NAMES OpenMM
        PATHS /usr/local/openmm/lib/
        DOC "openmm library"
    )
    find_library(
        OPENMM_KBPLUGIN
        NAMES KBPlugin
        PATHS /usr/local/openmm/lib/
        DOC "KBPlugin Plugin"
    )
    find_library(
        OPENMM_LIBRARY_CUDA
        NAMES OpenMMCUDA
        PATHS /usr/local/openmm/plugins/
        DOC "openmm library cuda"
    )
endif(CYGWIN)

set(OPENMM_INCLUDE_DIR ${OPENMM_INCLUDE_DIR} CACHE STRING INTERNAL)
set(OPENMM_LIBRARY ${OPENMM_LIBRARY} ${OPENMM_LIBRARY_CUDA} ${OPENMM_KBPLUGIN} CACHE STRING INTERNAL)

IF(OPENMM_LIBRARY AND OPENMM_INCLUDE_DIR AND OPENMM_KBPLUGIN)

SET(OPENMM_FOUND 1)

ENDIF(OPENMM_LIBRARY AND OPENMM_INCLUDE_DIR AND OPENMM_KBPLUGIN)

# ==========================================
IF(NOT OPENMM_FOUND)
    MESSAGE(FATAL_ERROR "OPENMM required, please specify it's location.")
ENDIF(NOT OPENMM_FOUND)
