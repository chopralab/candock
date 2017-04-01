include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( pdbreader_SRC
     pdbreader.cpp
     pdbreader.hpp
)

add_library(pdbreader OBJECT ${pdbreader_SRC})
