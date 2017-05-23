#.rst:
# FindOPENMM
# --------
#
# Find the native OPENMM includes and libraries.
#
# The GNU Scientific Library (OPENMM) is a numerical library for C and C++
# programmers. It is free software under the GNU General Public
# License.
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If OPENMM is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
#  OPENMM::openmm      - The main OPENMM library.
#  OPENMM::kbforce     - The KBFORCE support library used by OPENMM.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  OPENMM_FOUND          - True if OPENMM found on the local system
#  OPENMM_INCLUDE_DIRS   - Location of OPENMM header files.
#  OPENMM_LIBRARIES      - The OPENMM libraries.
#
# Hints
# ^^^^^
#
# Set ``OPENMM_ROOT_DIR`` to a directory that contains a OPENMM installation.
#
# This script expects to find libraries at ``$OPENMM_ROOT_DIR/lib`` and the OPENMM
# headers at ``$OPENMM_ROOT_DIR/include/``.  The library directory may
# optionally provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type
# of OPENMM installation discovered.  These variables may optionally be set to
# help this module find the correct files::
#
#  OPENMM_KBFORCE_LIBRARY       - Location of the OPENMM KBFORCE library.
#  OPENMM_KBFORCE_LIBRARY_DEBUG - Location of the debug OPENMM KBFORCE library (if any).
#  OPENMM_LIBRARY             - Location of the OPENMM library.
#  OPENMM_LIBRARY_DEBUG       - Location of the debug OPENMM library (if any).
#

#=============================================================================
# Copyright 2017 Jonathan Fine <finej@purdue.edu>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``OPENMM_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{OPENMM_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{OPENMM_ROOT_DIR}" OPENMM_ROOT_DIR )
  set( OPENMM_ROOT_DIR "${OPENMM_ROOT_DIR}" CACHE PATH "Prefix for OPENMM installation." )
endif()

#=============================================================================
# Set OPENMM_INCLUDE_DIRS and OPENMM_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $OPENMM_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( OPENMM_INCLUDE_DIR
  NAMES OpenMM.h
  HINTS ${OPENMM_ROOT_DIR}/include
)
find_library( OPENMM_LIBRARY
  NAMES OpenMM
  HINTS ${OPENMM_ROOT_DIR}/lib
  PATH_SUFFIXES Release Debug
)
find_library( OPENMM_KBFORCE_LIBRARY
  NAMES KBPlugin
  HINTS ${OPENMM_ROOT_DIR}/lib
  PATH_SUFFIXES Release Debug
)
# Do we also have debug versions?
find_library( OPENMM_LIBRARY_DEBUG
  NAMES OpenMM
  HINTS ${OPENMM_ROOT_DIR}/lib
  PATH_SUFFIXES Debug
)
find_library( OPENMM_KBFORCE_LIBRARY_DEBUG
  NAMES KBPlugin
  HINTS ${OPENMM_ROOT_DIR}/lib
  PATH_SUFFIXES Debug
)
set( OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIR} )
set( OPENMM_LIBRARIES ${OPENMM_LIBRARY} ${OPENMM_KBFORCE_LIBRARY} )

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set OPENMM_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( OPENMM
  FOUND_VAR
    OPENMM_FOUND
  REQUIRED_VARS
    OPENMM_INCLUDE_DIR
    OPENMM_LIBRARY
    OPENMM_KBFORCE_LIBRARY
    )

mark_as_advanced( OPENMM_ROOT_DIR OPENMM_LIBRARY OPENMM_INCLUDE_DIR
  OPENMM_KBFORCE_LIBRARY OPENMM_LIBRARY_DEBUG OPENMM_KBFORCE_LIBRARY_DEBUG
  OPENMM_USE_PKGCONFIG OPENMM_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" OPENMM_LIBRARY_DLL       "${OPENMM_LIBRARY}" )
  string( REPLACE ".lib" ".dll" OPENMM_KBFORCE_LIBRARY_DLL "${OPENMM_KBFORCE_LIBRARY}" )
  string( REPLACE ".lib" ".dll" OPENMM_LIBRARY_DEBUG_DLL "${OPENMM_LIBRARY_DEBUG}" )
  string( REPLACE ".lib" ".dll" OPENMM_KBFORCE_LIBRARY_DEBUG_DLL "${OPENMM_KBFORCE_LIBRARY_DEBUG}" )
endif()

if( OPENMM_FOUND AND NOT TARGET OPENMM::openmm )
  if( EXISTS "${OPENMM_LIBRARY_DLL}" AND EXISTS "${OPENMM_KBFORCE_LIBRARY_DLL}")

    # Windows systems with dll libraries.
    add_library( OPENMM::openmm      SHARED IMPORTED )
    add_library( OPENMM::kbforce SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( OPENMM::kbforce PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${OPENMM_KBFORCE_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${OPENMM_KBFORCE_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${OPENMM_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
    set_target_properties( OPENMM::openmm PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${OPENMM_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${OPENMM_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${OPENMM_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      INTERFACE_LINK_LIBRARIES          OPENMM::kbforce )

    # If we have both Debug and Release libraries
    if( EXISTS "${OPENMM_LIBRARY_DEBUG_DLL}" AND EXISTS "${OPENMM_KBFORCE_LIBRARY_DEBUG_DLL}")
      set_property( TARGET OPENMM::kbforce APPEND PROPERTY IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( OPENMM::kbforce PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${OPENMM_KBFORCE_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${OPENMM_KBFORCE_LIBRARY_DEBUG}" )
      set_property( TARGET OPENMM::openmm APPEND PROPERTY IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( OPENMM::openmm PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${OPENMM_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${OPENMM_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create
    # the imported library targets.
    add_library( OPENMM::openmm      UNKNOWN IMPORTED )
    add_library( OPENMM::kbforce UNKNOWN IMPORTED )
    set_target_properties( OPENMM::kbforce PROPERTIES
      IMPORTED_LOCATION                 "${OPENMM_KBFORCE_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${OPENMM_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
    set_target_properties( OPENMM::openmm PROPERTIES
      IMPORTED_LOCATION                 "${OPENMM_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${OPENMM_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      INTERFACE_LINK_LIBRARIES          OPENMM::kbforce )
  endif()
endif()
