include_directories($CMAKE_CURRENT_SOURCE_DIR)

add_subdirectory(centro)
add_subdirectory(cluster)
add_subdirectory(docker)
add_subdirectory(design)
add_subdirectory(fragmenter)
add_subdirectory(geom3d)
add_subdirectory(graph)
add_subdirectory(helper)
add_subdirectory(kabsch)
add_subdirectory(kbforce)
add_subdirectory(ligands)
add_subdirectory(linker)
add_subdirectory(modeler)
add_subdirectory(pdbreader)
add_subdirectory(probis)
add_subdirectory(score)
add_subdirectory(program)

if (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})

add_library(
    candock SHARED 
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

install(
    TARGETS
        candock
    LIBRARY DESTINATION
        lib
)

endif (NOT ${CANDOCK_BUILD_STATIC_EXECUTABLE})