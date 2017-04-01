include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( geom3d_SRC
        coordinate.cpp
        coordinate.hpp
        geom3d.cpp
        geom3d.hpp
        matrix.cpp
        matrix.hpp
        pca.cpp
        pca.hpp
        quaternion.hpp
)

add_library(geom3d OBJECT ${geom3d_SRC})
