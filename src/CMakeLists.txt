include_directories($CMAKE_CURRENT_SOURCE_DIR)

include_directories( ../foreign/tclap/include/ )
include_directories( ../lib/modeler/ )

set( candock_libraries 
        program
        docker
        score
        centro
        kabsch
        ligands
        cluster   # Ugly GCC hack to fix interdependency issue
        linker
        cluster
        graph
        modeler
        pdbreader
        probis
        fragmenter
        kbforce
        geom3d
        helper
   )

set( additional_libraries 
        pthread   # Required for GCC
        dl
   )

set( all_needed_libraries
     ${candock_libraries}
     ${Boost_LIBRARIES}
     ${OPENMM_LIBRARY}
     ${GSL_LIBRARIES}
     ${additional_libraries}
   )

add_executable( candock candock.cpp )
add_executable( genpot  genpot.cpp )
add_executable( test_minimization test_minimization.cpp )
add_executable( test_fragmenting test_fragmenting.cpp )
add_executable( test_bsite test_bsite.cpp )
add_executable( test_clq test_clq.cpp )
add_executable( test_link test_link.cpp )
add_executable( test_singlepoint test_singlepoint.cpp )
add_executable( test_cofactor test_cofactor.cpp )
add_executable( candock_oo candock_oo.cpp )

target_link_libraries( candock ${all_needed_libraries} )
target_link_libraries( genpot ${all_needed_libraries} )
target_link_libraries( test_minimization ${all_needed_libraries} )
target_link_libraries( test_fragmenting ${all_needed_libraries} )
target_link_libraries( test_bsite ${all_needed_libraries} )
target_link_libraries( test_clq ${all_needed_libraries} )
target_link_libraries( test_link ${all_needed_libraries} )
target_link_libraries( test_singlepoint ${all_needed_libraries} )
target_link_libraries( test_cofactor ${all_needed_libraries} )
target_link_libraries( candock_oo ${all_needed_libraries} design)
