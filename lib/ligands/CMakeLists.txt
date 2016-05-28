include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( ligands_SRC common.cpp
        				 common.hpp
				         gbio.cpp
				         genbio.cpp
				         genbio.hpp
				         genclus.cpp
				         genclus.hpp
				         genlig.cpp
				         genlig.hpp
				         glig.cpp
				         jsonreader.cpp
				         jsonreader.hpp
				         nosqlreader.cpp
				         nosqlreader.hpp
				         opts_genbio.hpp
				         opts_genclus.hpp
				         opts_genlig.hpp )

include_directories(../tclap-1.2.1/include )

add_library(ligands STATIC ${ligands_SRC})
