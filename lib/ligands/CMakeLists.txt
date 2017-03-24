include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( ligands_SRC
        common.cpp
        common.hpp
        jsonreader.cpp
        jsonreader.hpp
        genbio.cpp
        genbio.hpp
        genclus.cpp
        genclus.hpp
        genlig.cpp
        genlig.hpp
        nosqlreader.cpp
        nosqlreader.hpp
)



add_library(ligands OBJECT ${ligands_SRC})
