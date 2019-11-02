include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( score_SRC
        interpolation.cpp
        interpolation.hpp
        score.cpp
        score.hpp
        xscore.hpp
        xscore.cpp
)

add_library(score OBJECT ${score_SRC})
add_dependencies(score gsl)

