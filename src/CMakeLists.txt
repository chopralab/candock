include_directories($CMAKE_CURRENT_SOURCE_DIR)

include_directories( ../foreign/tclap-1.2.1/include )
include_directories( ../lib/modeler/ )

set( candock_libraries 
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

add_executable( candock candock.cpp common.cpp common.hpp opts_candock.hpp )
add_executable( genpot  genpot.cpp )
add_executable( test_minimization test_minimization.cpp common.cpp )
add_executable( test_fragmenting test_fragmenting.cpp common.cpp )
add_executable( test_bsite test_bsite.cpp common.cpp )
add_executable( test_clq test_clq.cpp common.cpp )
add_executable( test_link test_link.cpp common.cpp )
add_executable( test_singlepoint test_singlepoint.cpp common.cpp )
add_executable( test_cofactor test_cofactor.cpp common.cpp )

target_link_libraries( candock
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${OPENMM_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( genpot
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_minimization
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${OPENMM_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_fragmenting
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_bsite
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_clq
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_link
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${OPENMM_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_singlepoint
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( test_cofactor
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${OPENMM_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )
