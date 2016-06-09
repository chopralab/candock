# Find OpenMM
#
# Find the OpenMM includes and library
# 
# if you nee to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  OpenMM_INCLUDE_DIRS, where to find header, etc.
#  OpenMM_LIBRARIES, the libraries needed to use OpenMM.
#  OpenMM_FOUND, If false, do not try to use OpenMM.

# only look in default directories
find_path(
	OpenMM_INCLUDE_DIR 
	NAMES OpenMM.h
	DOC "OpenMM include dir"
)

find_library(
	OpenMM_LIBRARY
	NAMES OpenMM
	DOC "OpenMM library"
)

set(OpenMM_INCLUDE_DIRS ${OpenMM_INCLUDE_DIR})
set(OpenMM_LIBRARIES ${OpenMM_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set OpenMM_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMM DEFAULT_MSG
	OpenMM_INCLUDE_DIR OpenMM_LIBRARY)
mark_as_advanced (OpenMM_INCLUDE_DIR OpenMM_LIBRARY)
