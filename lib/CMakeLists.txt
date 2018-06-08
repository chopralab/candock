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

    if( ${CANDOCK_DRM} )
        set( candock_libraries_to_link ${candock_libraries_to_link} ${OPENSSL_LIBRARIES})
    endif( ${CANDOCK_DRM} )

    target_link_libraries ( candock ${candock_libraries_to_link} )

    if (NOT WIN32)
        install(
            TARGETS
                candock
            LIBRARY DESTINATION
                ${CANDOCK_VERSION}/lib
        )
    else()
        install(
            TARGETS
               candock
               ${GSL_LIBRARIES}
               ${OPENMM_LIBRARIES}
               ${Boost_LIBRARIES}   
            RUNTIME DESTINATION
                ${CANDOCK_VERSION}/bin/
        )
    endif()
else()
    add_library(
        candock STATIC
            candockexport.hpp
            ${CANDOCK_LIB_SOURCES}
    )
endif (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})

if (${USING_INTERNAL_GSL})
    add_dependencies(candock GSL)
endif()

if (${USING_INTERNAL_OPENMM})
    add_dependencies(candock OpenMM_Build)
endif()
