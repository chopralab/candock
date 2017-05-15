include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( drm_SRC
        drm.cpp
        drm.hpp
)

add_library(drm OBJECT ${drm_SRC})
