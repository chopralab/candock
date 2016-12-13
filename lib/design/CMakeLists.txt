include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( design_SRC design.cpp
				design.hpp )

add_library(design STATIC ${design_SRC})
