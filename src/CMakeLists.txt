
if(${CANDOCK_BUILD_STATIC_EXECUTABLE})

    if (WIN32)
            add_definitions( /DCANDOCK_USE_STATIC_LIBRARIES)
    endif(WIN32)

    set( CANDOCK_SRC
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

set( program_names  check_options
                    check_inputpdb
                    find_centroids
                    prep_fragments
                    make_objective
                    make_bsitegrid
                    dock_fragments
                    link_fragments
                    design_ligands
                    convertto_mol2
                    calc_rmsdcryst
                    calc_rmsdtable
                    calc_potential
                    calc_posescore
                    calc_posexscore
                    calc_allscores
                    minimize_poses
                    dynamics_run
)

foreach( candock_program ${program_names} )

    add_executable( cd_${candock_program} "cd_${candock_program}.cpp"  ${CANDOCK_SRC})

    target_link_libraries( cd_${candock_program} ${all_needed_libraries} )

    install(
        TARGETS
            cd_${candock_program}

        RUNTIME DESTINATION
        ${CANDOCK_VERSION}/bin
    )

    if( NOT WIN32 )
        install(CODE "message( STATUS \"Linking ${candock_program}.sh\")
            execute_process(
                  WORKING_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/modules
                  COMMAND ${CMAKE_COMMAND} -E create_symlink module.sh ${candock_program}.sh)")
    endif( NOT WIN32 )
endforeach( candock_program ${program_names} )

if( ${CANDOCK_DRM} )
    add_executable( drm_key_generator drm_key_generator.cpp )
    target_link_libraries( drm_key_generator ${OPENSSL_LIBRARIES} )
    install( TARGETS drm_key_generator RUNTIME DESTINATION ${CANDOCK_VERSION}/bin )
endif( ${CANDOCK_DRM} )
