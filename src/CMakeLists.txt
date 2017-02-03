include_directories($CMAKE_CURRENT_SOURCE_DIR)

include_directories( ../foreign/tclap/include/ )
include_directories( ../lib/modeler/ )

set( candock_libraries 
        program
        design
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
     /usr/local/openmm/lib/plugins/libOpenMMCUDA.so
   )

# Original, to be deprecated programs
add_executable( candock candock.cpp )
add_executable( genpot  genpot.cpp )
add_executable( test_minimization test_minimization.cpp )
add_executable( test_fragmenting test_fragmenting.cpp )
add_executable( test_bsite test_bsite.cpp )
add_executable( test_clq test_clq.cpp )
add_executable( test_link test_link.cpp )
add_executable( test_singlepoint test_singlepoint.cpp )
add_executable( test_cofactor test_cofactor.cpp )

# New programs with new features that change functionality
add_executable( cd_check_options  cd_check_options.cpp )
add_executable( cd_bsite_probis   cd_bsite_probis.cpp )
add_executable( cd_prep_fragments cd_prep_fragments.cpp )
add_executable( cd_dock_fragments cd_dock_fragments.cpp )
add_executable( cd_link_fragments cd_link_fragments.cpp )
add_executable( cd_design_ligands cd_design_ligands.cpp )

if(CUDA_FOUND)
#CUDA programs
cuda_add_executable( cuda_link_fragments cd_link_fragments.cpp )
endif(CUDA_FOUND)

target_link_libraries( candock ${all_needed_libraries} )
target_link_libraries( genpot ${all_needed_libraries} )
target_link_libraries( test_minimization ${all_needed_libraries} )
target_link_libraries( test_fragmenting ${all_needed_libraries} )
target_link_libraries( test_bsite ${all_needed_libraries} )
target_link_libraries( test_clq ${all_needed_libraries} )
target_link_libraries( test_link ${all_needed_libraries} )
target_link_libraries( test_singlepoint ${all_needed_libraries} )
target_link_libraries( test_cofactor ${all_needed_libraries} )

target_link_libraries( cd_check_options  ${all_needed_libraries} )
target_link_libraries( cd_bsite_probis   ${all_needed_libraries} )
target_link_libraries( cd_prep_fragments ${all_needed_libraries} )
target_link_libraries( cd_dock_fragments ${all_needed_libraries} )
target_link_libraries( cd_link_fragments ${all_needed_libraries} )
target_link_libraries( cd_design_ligands ${all_needed_libraries} )

if(CUDA_FOUND)
target_link_libraries( cuda_link_fragments ${all_needed_libraries} )
endif(CUDA_FOUND)

install( TARGETS
    candock genpot
    test_minimization
    test_fragmenting
    test_bsite
    test_clq
    test_link
    test_singlepoint
    test_cofactor
    
    cd_check_options
    cd_bsite_probis
    cd_prep_fragments
    cd_dock_fragments
    cd_link_fragments
    cd_design_ligands
    
    RUNTIME DESTINATION v${CANDOCK_MAJOR_VERSION}.${CANDOCK_MINOR_VERSION}.${CANDOCK_TWEAK_VERSION}
)
