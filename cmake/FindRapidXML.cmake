# This script locates the RapidXML library
# ------------------------------------
#
# usage:
# find_package(RAPIDXML ...)
#
# searches in RAPIDXML_ROOT and usual locations
#
# Sets RAPIDXML_DIR

find_path(RAPIDXML_DIR 
          NAMES rapidxml/rapidxml.hpp
          DOC "RapidXML include dir"
          HINTS ${RAPIDXML_ROOT}
)

if (NOT RAPIDXML_DIR)
    if(RapidXML_FIND_REQUIRED) #prefix is filename, case matters
        message(FATAL_ERROR "Could not find RapidXML!")
    elseif(NOT RapidXML_FIND_QUIETLY)
        message("Could not find RapidXML!")
    endif(RapidXML_FIND_REQUIRED)
endif(NOT RAPIDXML_DIR)

if (NOT RapidXML_FIND_QUIETLY)
    message("Found RapidXML: ${RAPIDXML_DIR}")
endif (NOT RapidXML_FIND_QUIETLY)
