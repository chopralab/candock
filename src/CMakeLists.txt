
if(${CANDOCK_BUILD_STATIC_EXECUTABLE})

    if (WIN32)
            add_definitions( /DCANDOCK_USE_STATIC_LIBRARIES)
    endif(WIN32)

    set( all_needed_libraries
        candock
        ${OPENMM_LIBRARY}
        ${Boost_LIBRARIES}
        ${GSL_LIBRARIES}
    )

    if ( ${CANDOCK_DRM} )
        set( all_needed_libraries ${all_needed_libraries} ${OPENSSL_LIBRARIES} )
    endif( ${CANDOCK_DRM} )

else (${CANDOCK_BUILD_STATIC_EXECUTABLE})
    
    set( all_needed_libraries
        candock
        ${STATCHEM_LIBRARY}
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
                    find_centroids
                    make_bsitegrid
                    dock_fragments
                    link_fragments
                    design_ligands
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
