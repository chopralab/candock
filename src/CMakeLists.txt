include_directories($CMAKE_CURRENT_SOURCE_DIR)

include_directories( ../foreign/tclap/include/ )
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

add_executable( candock candock.cpp common.cpp common.hpp opts_candock.cpp )
add_executable( genpot  genpot.cpp opts_candock.cpp )
add_executable( test_minimization test_minimization.cpp common.cpp opts_candock.cpp )
add_executable( test_fragmenting test_fragmenting.cpp common.cpp opts_candock.cpp )
add_executable( test_bsite test_bsite.cpp common.cpp opts_candock.cpp )
add_executable( test_clq test_clq.cpp common.cpp opts_candock.cpp )
add_executable( test_link test_link.cpp common.cpp opts_candock.cpp)
add_executable( test_singlepoint test_singlepoint.cpp common.cpp opts_candock.cpp)
add_executable( test_cofactor test_cofactor.cpp common.cpp opts_candock.cpp)

add_executable( candock_oo candock_oo.cpp
	common.hpp
	common.cpp
	opts_candock.hpp
	opts_candock.cpp
	programstep.hpp
	findcentroids.hpp
	findcentroids.cpp
	fragmentligands.hpp
	fragmentligands.cpp
	dockfragments.hpp
	dockfragments.cpp
	linkfragments.hpp
	linkfragments.cpp
	target.hpp
	target.cpp
)

target_link_libraries( candock
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
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
                         ${Boost_program_options_LIBRARY}
                         ${OPENMM_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                     )

target_link_libraries( candock_oo
                         ${candock_libraries}
                         ${Boost_LIBRARIES}
                         ${Boost_date_time_LIBRARY}
                         ${Boost_filesystem_LIBRARY}
                         ${Boost_regex_LIBRARY}
                         ${Boost_system_LIBRARY}
                         ${Boost_program_options_LIBRARY}
                         ${OPENMM_LIBRARY}
                         ${GSL_LIBRARIES}
                         ${additional_libraries}
                         design
                         score # another ugly hack???
                     )
