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
   )

# New programs with new features that change functionality
add_executable( cd_check_options  cd_check_options.cpp )
add_executable( cd_bsite_probis   cd_bsite_probis.cpp )
add_executable( cd_prep_fragments cd_prep_fragments.cpp )
add_executable( cd_make_objective cd_make_objective.cpp )
add_executable( cd_dock_fragments cd_dock_fragments.cpp )
add_executable( cd_link_fragments cd_link_fragments.cpp )
add_executable( cd_design_ligands cd_design_ligands.cpp )

if(CUDA_FOUND)
#CUDA programs
cuda_add_executable( cuda_link_fragments cd_link_fragments.cpp )
endif(CUDA_FOUND)

target_link_libraries( cd_check_options  ${all_needed_libraries} )
target_link_libraries( cd_bsite_probis   ${all_needed_libraries} )
target_link_libraries( cd_prep_fragments ${all_needed_libraries} )
target_link_libraries( cd_make_objective ${all_needed_libraries} )
target_link_libraries( cd_dock_fragments ${all_needed_libraries} )
target_link_libraries( cd_link_fragments ${all_needed_libraries} )
target_link_libraries( cd_design_ligands ${all_needed_libraries} )

if(CUDA_FOUND)
target_link_libraries( cuda_link_fragments ${all_needed_libraries} )
endif(CUDA_FOUND)

install( TARGETS    
    cd_check_options
    cd_bsite_probis
    cd_prep_fragments
    cd_make_objective
    cd_dock_fragments
    cd_link_fragments
    cd_design_ligands
    
    RUNTIME DESTINATION v${CANDOCK_MAJOR_VERSION}.${CANDOCK_MINOR_VERSION}.${CANDOCK_TWEAK_VERSION}
)
