include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( cluster_SRC greedy.cpp
                 greedy.hpp
                 optics.hpp )

add_library(cluster STATIC ${cluster_SRC})
