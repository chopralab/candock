include_directories($CMAKE_CURRENT_SOURCE_DIR)
include_directories(${OPENMM_KBPLUGIN_INCLUDE})

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

if (USING_INTERNAL_OPENMM)
    add_dependencies(modeler OpenMM_Build)
endif()