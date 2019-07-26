# An external project for OpenMM
set(openmm_source  "${CMAKE_CURRENT_BINARY_DIR}/openmm-src")
set(openmm_install "${CMAKE_BINARY_DIR}/stage")
set(openmm_version "7.2.2")
set(openmm_url "https://github.com/pandegroup/openmm/archive/${openmm_version}.tar.gz")
set(openmm_md5 "06c44b6703c6f1d6f2908d56600ad6b6")

ExternalProject_Add(openmm
  DOWNLOAD_DIR ${download_dir}
  SOURCE_DIR "${openmm_source}"
  INSTALL_DIR "${openmm_install}"
  URL ${openmm_url}
  URL_MD5 ${openmm_md5}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  CMAKE_CACHE_ARGS
    -DOPENMM_BUILD_C_AND_FORTRAN_WRAPPERS:Bool=off
    -DOPENMM_BUILD_PYTHON_WRAPPERS:Bool=off
    -DOPENMM_BUILD_EXAMPLES:Bool=off
    -DOPENMM_BUILD_TESTS:Bool=off
    -DOPENMM_BUILD_REFERENCE_TESTS:Bool=off
    -DOPENMM_BUILD_AMOEBA_PLUGIN:Bool=on
    -DOPENMM_BUILD_RPMD_PLUGIN:Bool=on
    -DOPENMM_BUILD_DRUDE_PLUGIN:Bool=on
    -DOPENMM_BUILD_PME_PLUGIN:Bool=off
    -DBUILD_TESTING:BOOL=OFF
)

set(OPENMM_LIBRARIES ${CMAKE_BINARY_DIR}/stage/lib/libOpenMM.so)

install(
    FILES
        ${CMAKE_BINARY_DIR}/stage/lib/libOpenMM.so
    DESTINATION
        lib/${CANDOCK_VERSION}
)


