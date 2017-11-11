add_executable( structure_mutator
                structure_mutator.cpp
            )

target_link_libraries( structure_mutator candock ${Boost_LIBRARIES} )

install(
    TARGETS
        structure_mutator

    RUNTIME DESTINATION
        bin/${CANDOCK_VERSION}/tools
)

if( NOT WIN32 )
    install(CODE "message( STATUS \"Linking gen_count_file.sh\")
        execute_process(
                WORKING_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/modules/${CANDOCK_VERSION}
                COMMAND ${CMAKE_COMMAND} -E create_symlink module.sh gen_count_file.sh)")
endif( NOT WIN32 )
