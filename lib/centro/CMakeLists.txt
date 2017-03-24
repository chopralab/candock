include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( centro_SRC
        centroids.cpp
        centroids.hpp
)

add_library(centro OBJECT ${centro_SRC})
