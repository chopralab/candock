include_directories($CMAKE_CURRENT_SOURCE_DIR)

if (WIN32)
    if   (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})
            add_definitions( /DCANDOCK_SHARED_LIBRARY)
    else (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})
            add_definitions( /DCANDOCK_STATIC_LIBRARY)
    endif(NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})
endif(WIN32)

file(GLOB_RECURSE CANDOCK_LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/**.cpp)

if (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})

    add_library(
        candock SHARED
            ${CANDOCK_LIB_SOURCES}
    )

    set( candock_libraries_to_link ${OPENMM_LIBRARY}
                                   ${Boost_LIBRARIES}
                                   ${GSL_LIBRARIES}
    )

    target_link_libraries ( candock ${candock_libraries_to_link} )
else()
    add_library(
        candock STATIC
            ${CANDOCK_LIB_SOURCES}
    )
endif()

install(
    TARGETS
        candock
    LIBRARY DESTINATION
        lib
    ARCHIVE DESTINATION
        lib
    RUNTIME DESTINATION
        bin
)
