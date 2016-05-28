include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( fragmenter_SRC fragmenter.cpp
		    fragmenter.hpp
		    unique.cpp
		    unique.hpp )

add_library(fragmenter STATIC ${fragmenter_SRC} )
