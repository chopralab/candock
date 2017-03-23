include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( score_SRC
        interpolation.cpp
        interpolation.hpp
        score.cpp
        score.hpp
)

add_library(score OBJECT ${score_SRC})
