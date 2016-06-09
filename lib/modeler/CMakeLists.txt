include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( modeler_SRC forcefield.cpp
        				 forcefield.hpp
			        	 modeler.cpp
				         modeler.hpp
				         systemtopology.cpp
				         systemtopology.hpp
				         topology.cpp
				         topology.hpp )

include_directories( ../modeler )

add_library(modeler STATIC ${modeler_SRC})
