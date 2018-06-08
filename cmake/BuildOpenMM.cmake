
set(OPENMM_GIT_HASH "8fe2b5a5c4c399cdebdf9e6950c2ead9a17b586f")

include(ExternalProject)

set(USING_INTERNAL_OPENMM on CACHE BOOL "Using Internal GSL")
mark_as_advanced(USING_INTERNAL_OPENMM)

ExternalProject_Add(OpenMM_Build
    URL https://github.com/pandegroup/openmm/archive/${OPENMM_GIT_HASH}.tar.gz
    URL_HASH SHA1=4c7c70e6113d31e11a9efe574e6626619e0de34b
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmm
    CMAKE_CACHE_ARGS
                     -DOPENMM_BUILD_C_AND_FORTRAN_WRAPPERS:Bool=off
                     -DOPENMM_BUILD_PYTHON_WRAPPERS:Bool=off
                     -DOPENMM_BUILD_EXAMPLES:Bool=off
                     -DOPENMM_BUILD_TESTS:Bool=off
                     -DOPENMM_BUILD_REFERENCE_TESTS:Bool=off
                     -DOPENMM_BUILD_AMOEBA_PLUGIN:Bool=off
                     -DOPENMM_BUILD_RPMD_PLUGIN:Bool=off
                     -DOPENMM_BUILD_DRUDE_PLUGIN:Bool=off
                     -DOPENMM_BUILD_PME_PLUGIN:Bool=off
                     -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/openmm_master
    CMAKE_ARGS -DBUILD_TESTING=off
)

add_library(OpenMM SHARED IMPORTED)
set_property(TARGET OpenMM PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/openmm_master/lib/${CMAKE_SHARED_LIBRARY_PREFIX}OpenMM${CMAKE_SHARED_LIBRARY_SUFFIX} )
add_dependencies(OpenMM OpenMM_Build)

set(OPENMM_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmm_master/include CACHE PATH "Internal OpenMM Includes")
set(OPENMM_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/openmm_master/lib/${CMAKE_SHARED_LIBRARY_PREFIX}OpenMM${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE STRING "Internal OpenMM Library" )

mark_as_advanced(OPENMM_LIBRARY)
mark_as_advanced(OPENMM_INCLUDE_DIR)

#install(FILES
# ${CMAKE_CURRENT_BINARY_DIR}/openmm_master/lib/${CMAKE_SHARED_LIBRARY_PREFIX}OpenMM${CMAKE_SHARED_LIBRARY_SUFFIX}
#DESTINATION
# ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/lib
#)

set(OPENMM_PLUGIN_DIR ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/lib/plugins CACHE STRING "OpenMM Plugin Location")
mark_as_advanced(OPENMM_PLUGIN_DIR)
