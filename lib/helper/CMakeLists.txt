include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( helper_SRC
        array1d.hpp
        array2d.cpp
        array2d.hpp
        array3d.hpp
        benchmark.cpp
        benchmark.hpp
        debug.hpp
        error.hpp
        grep.cpp
        grep.hpp
        help.cpp
        help.hpp
        inout.cpp
        inout.hpp
        options.cpp
        options.hpp
        path.cpp
        path.hpp
        renamerules.hpp
        smiles.cpp
        smiles.hpp
)

add_library(helper OBJECT ${helper_SRC})
