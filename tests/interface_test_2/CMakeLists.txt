if( NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE} )

    add_executable(test_interface_2 test_interface_2.c)

    target_link_libraries(test_interface_2 candock ${Boost_LIBRARIES}
        ${OPENMM_LIBRARY}
        ${GSL_LIBRARIES})
endif( NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE} )