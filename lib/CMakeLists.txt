include_directories($CMAKE_CURRENT_SOURCE_DIR)

if (WIN32)
    if   (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})
            add_definitions( /DCANDOCK_SHARED_LIBRARY)
    else (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})
            add_definitions( /DCANDOCK_STATIC_LIBRARY)
    endif(NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})

endif(WIN32)

add_subdirectory(centro)
add_subdirectory(cluster)
add_subdirectory(drm)
add_subdirectory(fileout)
add_subdirectory(docker)
add_subdirectory(design)
add_subdirectory(fragmenter)
add_subdirectory(geom3d)
add_subdirectory(graph)
add_subdirectory(helper)
add_subdirectory(interface)
add_subdirectory(kabsch)
add_subdirectory(ligands)
add_subdirectory(linker)
add_subdirectory(modeler)
add_subdirectory(molib)
add_subdirectory(parser)
add_subdirectory(probis)
add_subdirectory(score)
add_subdirectory(program)

if (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})

add_library(
    candock SHARED
        candockexport.hpp
        $<TARGET_OBJECTS:interface>
        $<TARGET_OBJECTS:fileout>
        $<TARGET_OBJECTS:program>
        $<TARGET_OBJECTS:design>
        $<TARGET_OBJECTS:docker>
        $<TARGET_OBJECTS:drm>
        $<TARGET_OBJECTS:score>
        $<TARGET_OBJECTS:centro>
        $<TARGET_OBJECTS:kabsch>
        $<TARGET_OBJECTS:ligands>
        $<TARGET_OBJECTS:linker>
        $<TARGET_OBJECTS:cluster>
        $<TARGET_OBJECTS:graph>
        $<TARGET_OBJECTS:modeler>
        $<TARGET_OBJECTS:parser>
        $<TARGET_OBJECTS:molib>
        $<TARGET_OBJECTS:probis>
        $<TARGET_OBJECTS:fragmenter>
        $<TARGET_OBJECTS:geom3d>
        $<TARGET_OBJECTS:helper>
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

else (NOT WIN32)

    install(
        TARGETS
            candock
        RUNTIME DESTINATION
        ${CANDOCK_VERSION}/bin/
    )

endif(NOT WIN32)

endif (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})

if(WIN32)

     install(
	     FILES
			${GSL_LIBRARIES}
			${OPENMM_LIBRARIES}
			${Boost_LIBRARIES}
	 	 RUNTIME DESTINATION
          ${CANDOCK_VERSION}/lib
	 )

endif(WIN32)
