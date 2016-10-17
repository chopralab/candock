include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( program_SRC common.cpp
			     common.hpp
			     cmdlnopts.cpp
			     cmdlnopts.hpp
			     programstep.hpp
			     findcentroids.cpp
			     findcentroids.hpp
			     fragmentligands.cpp
			     fragmentligands.hpp
			     dockfragments.cpp
			     dockfragments.hpp
			     linkfragments.cpp
			     linkfragments.hpp)

add_library(program STATIC ${program_SRC})
