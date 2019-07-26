include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( linker_SRC
        dockedconformation.cpp
        dockedconformation.hpp
        genericlinker.cpp
        iterativelinker.cpp
        linker.cpp
        linker.hpp
        partial.cpp
        partial.hpp
        poses.cpp
        poses.hpp
        seed.cpp
        seed.hpp
        segment.cpp
        segment.hpp
        state.cpp
        state.hpp
        staticlinker.cpp
)

include_directories( ../modeler )

add_library(linker OBJECT ${linker_SRC})
add_dependencies(linker kbforce gsl openmm)

