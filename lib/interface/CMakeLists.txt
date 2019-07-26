include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( interface_SRC
        interface.cpp
        interface.hpp)

add_library(interface OBJECT ${interface_SRC})
add_dependencies(interface kbforce gsl openmm)

