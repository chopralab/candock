include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( pdbreader_SRC
     parser.cpp
     parser.hpp
     pdbreader.cpp
     pdbreader.hpp
)

add_library(pdbreader OBJECT ${pdbreader_SRC})
