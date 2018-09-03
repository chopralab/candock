add_executable( structure_mutator
                structure_mutator.cpp
            )

target_link_libraries( structure_mutator
    candock
    ${STATCHEM_LIBRARY}
    ${Boost_LIBRARIES}
    ${OPENMM_LIBRARY}
    ${GSL_LIBRARIES}
)

install(
    TARGETS
        structure_mutator

    RUNTIME DESTINATION
        ${CANDOCK_VERSION}/bin
)

if( NOT WIN32 )
    install(CODE "message( STATUS \"Linking structure_mutator.sh\")
        execute_process(
                WORKING_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/modules
                COMMAND ${CMAKE_COMMAND} -E create_symlink module.sh structure_mutator.sh)")
endif( NOT WIN32 )
