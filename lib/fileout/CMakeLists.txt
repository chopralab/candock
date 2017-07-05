include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( fileout_SRC
    fileout.hpp
    complex.cpp
    print_mol2.cpp
)

add_library(fileout OBJECT ${fileout_SRC})
