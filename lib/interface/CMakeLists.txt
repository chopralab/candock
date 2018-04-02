include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( interface_SRC
        interface.cpp
        interface.hpp)

add_library(interface OBJECT ${interface_SRC})

if (USING_INTERNAL_OPENMM)
    add_dependencies(interface OpenMM_Build)
endif()