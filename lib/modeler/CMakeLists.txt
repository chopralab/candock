include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( modeler_SRC
        forcefield.cpp
        forcefield.hpp
        modeler.cpp
        modeler.hpp
        systemtopology.cpp
        systemtopology.hpp
        topology.cpp
        topology.hpp
)

add_library(modeler OBJECT ${modeler_SRC})
add_dependencies(modeler kbforce gsl openmm)

