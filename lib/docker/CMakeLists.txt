include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( docker_SRC
        conformations.cpp
        conformations.hpp
        dock.cpp
        dock.hpp
        gpoints.cpp
        gpoints.hpp
)

add_library(docker OBJECT ${docker_SRC})
