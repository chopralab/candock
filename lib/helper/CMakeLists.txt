include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( helper_SRC array1d.hpp
		        		array2d.cpp
				        array2d.hpp
				        array3d.hpp
				        benchmark.cpp
				        benchmark.hpp
				        debug.hpp
				        error.hpp
				        help.cpp
				        help.hpp
				        inout.cpp
				        inout.hpp
				        path.cpp
				        path.hpp )

add_library(helper STATIC ${helper_SRC})
