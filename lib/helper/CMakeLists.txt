include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( helper_SRC
        array1d.hpp
        array2d.hpp
        array3d.hpp
        benchmark.hpp
        debug.hpp
        error.hpp
        grep.cpp
        grep.hpp
        help.cpp
        help.hpp
        inout.cpp
        inout.hpp
        path.cpp
        path.hpp
        renamerules.hpp
        smiles.cpp
        smiles.hpp
)

add_library(helper OBJECT ${helper_SRC})
