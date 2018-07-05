if (WIN32)
    if   (${CANDOCK_BUILD_STATIC_EXECUTABLE})
            add_definitions( /DCANDOCK_USE_STATIC_LIBRARIES)
    endif(${CANDOCK_BUILD_STATIC_EXECUTABLE})
endif(WIN32)

add_subdirectory( interface_test )
add_subdirectory( interface_test_2 )
add_subdirectory( interface_test_3 )
