add_executable( cd_gen_count_file
                cd_gen_count_file.cpp
                atominfotypes.hpp
                molecularbondextractor.hpp
                molecularbondextractor.cpp
                molecularbondbinner.hpp
                molecularbondbinner.cpp
            )

target_link_libraries( cd_gen_count_file candock ${Boost_LIBRARIES} )

install(
    TARGETS
        cd_gen_count_file

    RUNTIME DESTINATION
        ${CANDOCK_VERSION}/bin
)

if( NOT WIN32 )
    install(CODE "message( STATUS \"Linking gen_count_file.sh\")
        execute_process(
                WORKING_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/modules
                COMMAND ${CMAKE_COMMAND} -E create_symlink module.sh gen_count_file.sh)")
endif( NOT WIN32 )
