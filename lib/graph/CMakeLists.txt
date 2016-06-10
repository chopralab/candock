include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( graph_SRC graph.hpp
	      		   mcqd.cpp
			         mcqd.hpp
			         mnts.cpp
			         mnts.hpp
			         ullsubstate.hpp
			         vertex.hpp )

add_library(graph STATIC ${graph_SRC})

