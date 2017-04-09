include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( interface_SRC
        interface.cpp
        interface.hpp)

add_library(interface OBJECT ${interface_SRC})

add_executable(test_interface test_interface.c)

target_link_libraries(test_interface candock ${Boost_LIBRARIES}
        ${OPENMM_LIBRARY}
        ${GSL_LIBRARIES})
