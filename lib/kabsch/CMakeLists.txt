include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( kabsch_SRC kabsch.cpp
		        		kabsch.hpp )

add_library(kabsch STATIC ${kabsch_SRC})
