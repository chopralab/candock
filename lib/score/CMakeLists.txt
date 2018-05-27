include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( score_SRC
        kbff.hpp
        kbff.cpp
        interpolation.cpp
        interpolation.hpp
        powerfit.cpp
        powerfit.hpp
        score.cpp
        score.hpp
        xscore.hpp
        xscore.cpp
)

add_library(score OBJECT ${score_SRC})

if (USING_INTERNAL_GSL)
    add_dependencies(score GSL)
endif()
