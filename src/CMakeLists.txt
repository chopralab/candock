
if(${CANDOCK_BUILD_STATIC_EXECUTABLE})

    if (WIN32)
            add_definitions( /DCANDOCK_USE_STATIC_LIBRARIES)
    endif(WIN32)

    set( CANDOCK_SRC
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
    
    set( all_needed_libraries
        ${OPENMM_LIBRARIES}
        ${OPENMM_KBPLUGIN}
        ${Boost_LIBRARIES}
        ${GSL_LIBRARIES}
    )

    if ( ${CANDOCK_DRM} )
        set( all_needed_libraries ${all_needed_libraries} ${OPENSSL_LIBRARIES} )
    endif( ${CANDOCK_DRM} )

else (${CANDOCK_BUILD_STATIC_EXECUTABLE})

    set( CANDOCK_SRC )
    
    set( all_needed_libraries
        candock
    )

endif (${CANDOCK_BUILD_STATIC_EXECUTABLE})

if (NOT WIN32)

    set ( all_needed_libraries
        ${all_needed_libraries}
        pthread
        dl
    )

endif(NOT WIN32)

set( program_names cd_check_options
                    cd_find_centroids
                    cd_prep_fragments
                    cd_make_objective
                    cd_dock_fragments
                    cd_link_fragments
                    cd_cacrystal_rmsd
                    cd_calculate_rmsd
                    cd_calc_potential
                    cd_rescore_docked
                    cd_minimize_force
                    cd_minimize_poses
)

# New programs with new features that change functionality

foreach( candock_program ${program_names} )

    add_executable( ${candock_program} "${candock_program}.cpp"  ${CANDOCK_SRC})

    target_link_libraries( ${candock_program}  ${all_needed_libraries} )

    set_target_properties(
        ${candock_program}
        PROPERTIES
        INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib/${CANDOCK_VERSION}:${CMAKE_INSTALL_PREFIX}/lib/ext
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )

    install(
        TARGETS
            ${candock_program}

        RUNTIME DESTINATION
            bin/${CANDOCK_VERSION}
    )

endforeach( candock_program ${program_names} )

if( ${CANDOCK_DRM} )
    add_executable( drm_key_generator drm_key_generator.cpp )
    target_link_libraries( drm_key_generator ${OPENSSL_LIBRARIES} )
    install( TARGETS drm_key_generator RUNTIME DESTINATION bin/${CANDOCK_VERSION} )
endif( ${CANDOCK_DRM} )
