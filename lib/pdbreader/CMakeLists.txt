include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( pdbreader_SRC
     parser.cpp
     parser.hpp
     fileparser.cpp
     fileparser.hpp
     pdbparser.cpp
     mol2parser.cpp
)

add_library(pdbreader OBJECT ${pdbreader_SRC})
