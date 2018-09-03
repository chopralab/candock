# Find the Statchem library
#
# Defines:
#
#  STATCHEM_FOUND        - system has STATCHEM
#  STATCHEM_INCLUDE_DIRS - the STATCHEM include directories
#  STATCHEM_LIBRARY      - The STATCHEM library
#
find_path(STATCHEM_INCLUDE_DIR statchem/statchemexport.hpp)
find_library(STATCHEM_LIBRARY NAMES statchem)

set(STATCHEM_INCLUDE_DIRS "${STATCHEM_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(STATCHEM DEFAULT_MSG STATCHEM_INCLUDE_DIR STATCHEM_LIBRARY)

mark_as_advanced(STATCHEM_INCLUDE_DIR STATCHEM_LIBRARY)
