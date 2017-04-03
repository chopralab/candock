include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( parser_SRC
     parser.cpp
     parser.hpp
     fileparser.cpp
     fileparser.hpp
     pdbparser.cpp
     mol2parser.cpp
)

add_library(parser OBJECT ${parser_SRC})
