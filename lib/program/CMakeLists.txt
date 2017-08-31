include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( program_SRC
        candockexport.hpp
        cmdlnopts.cpp
        cmdlnopts.hpp
        programstep.hpp
        findcentroids.cpp
        findcentroids.hpp
        fragmentligands.cpp
        fragmentligands.hpp
        dockfragments.cpp
        dockfragments.hpp
        options.cpp
        options.hpp
        linkfragments.cpp
        linkfragments.hpp
        target.cpp
        target.hpp
        targetgroup.cpp
        targetgroup.hpp
  )

add_library(program OBJECT ${program_SRC})
