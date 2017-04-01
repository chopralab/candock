
if(${CANDOCK_BUILD_STATIC_EXECUTABLE})

    set( CANDOCK_SRC
        $<TARGET_OBJECTS:program>
        $<TARGET_OBJECTS:design>
        $<TARGET_OBJECTS:docker>
        $<TARGET_OBJECTS:score>
        $<TARGET_OBJECTS:centro>
        $<TARGET_OBJECTS:kabsch>
        $<TARGET_OBJECTS:ligands>
        $<TARGET_OBJECTS:linker>
        $<TARGET_OBJECTS:cluster>
        $<TARGET_OBJECTS:graph>
        $<TARGET_OBJECTS:modeler>
        $<TARGET_OBJECTS:pdbreader>
        $<TARGET_OBJECTS:probis>
        $<TARGET_OBJECTS:fragmenter>
        $<TARGET_OBJECTS:kbforce>
        $<TARGET_OBJECTS:geom3d>
        $<TARGET_OBJECTS:helper>
    )
    
    set( all_needed_libraries
        ${Boost_LIBRARIES}
        ${OPENMM_LIBRARY}
        ${GSL_LIBRARIES}
    )

else (${CANDOCK_BUILD_STATIC_EXECUTABLE})

    set( CANDOCK_SRC )
    
    set( all_needed_libraries
        candock
        ${Boost_LIBRARIES}
        ${OPENMM_LIBRARY}
        ${GSL_LIBRARIES}
    )

endif (${CANDOCK_BUILD_STATIC_EXECUTABLE})

if (NOT WIN32)

set ( all_needed_libraries
    ${all_needed_libraries}
    pthread
    dl
)

endif(NOT WIN32)

# New programs with new features that change functionality
add_executable( cd_check_options  cd_check_options.cpp  ${CANDOCK_SRC})
add_executable( cd_bsite_probis   cd_bsite_probis.cpp   ${CANDOCK_SRC})
add_executable( cd_prep_fragments cd_prep_fragments.cpp ${CANDOCK_SRC})
add_executable( cd_make_objective cd_make_objective.cpp ${CANDOCK_SRC})
add_executable( cd_dock_fragments cd_dock_fragments.cpp ${CANDOCK_SRC})
add_executable( cd_link_fragments cd_link_fragments.cpp ${CANDOCK_SRC})
add_executable( cd_rescore_docked cd_rescore_docked.cpp ${CANDOCK_SRC})
add_executable( cd_design_ligands cd_design_ligands.cpp ${CANDOCK_SRC})

target_link_libraries( cd_check_options  ${all_needed_libraries} )
target_link_libraries( cd_bsite_probis   ${all_needed_libraries} )
target_link_libraries( cd_prep_fragments ${all_needed_libraries} )
target_link_libraries( cd_make_objective ${all_needed_libraries} )
target_link_libraries( cd_dock_fragments ${all_needed_libraries} )
target_link_libraries( cd_link_fragments ${all_needed_libraries} )
target_link_libraries( cd_rescore_docked ${all_needed_libraries} )
target_link_libraries( cd_design_ligands ${all_needed_libraries} )

set_target_properties(
        cd_check_options
        cd_bsite_probis
        cd_prep_fragments
        cd_make_objective
        cd_dock_fragments
        cd_link_fragments
        cd_design_ligands
		cd_rescore_docked
    PROPERTIES
    INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib
    INSTALL_RPATH_USE_LINK_PATH TRUE
)


if(CUDA_FOUND)
    cuda_add_executable( cuda_link_fragments cd_link_fragments.cpp )
    target_link_libraries( cuda_link_fragments ${all_needed_libraries} )
    set_target_properties(
        cuda_link_fragments
        PROPERTIES
        INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
endif(CUDA_FOUND)

install(
    TARGETS    
        cd_check_options
        cd_bsite_probis
        cd_prep_fragments
        cd_make_objective
        cd_dock_fragments
        cd_link_fragments
        cd_design_ligands
		cd_rescore_docked

    RUNTIME DESTINATION
        bin
)
