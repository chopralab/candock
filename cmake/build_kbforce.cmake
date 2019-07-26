# An external project for OpenMM
set(kbforce_install "${CMAKE_BINARY_DIR}/stage")

ExternalProject_Add(kbforce
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/kbforce"
  INSTALL_DIR "${kbforce_install}"
  CMAKE_ARGS
     -DOPENMM_ROOT_DIR:PATH=${CMAKE_BINARY_DIR}/stage
     -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  #  -DCMAKE_BUILD_TYPE:STRING=Release
  #  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  #  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  DEPENDS
     openmm
)

set(OPENMM_KBPLUGIN ${CMAKE_BINARY_DIR}/stage/lib/libKBPlugin.so)

install(
    FILES
        ${CMAKE_BINARY_DIR}/stage/lib/libKBPlugin.so
    DESTINATION
        lib/${CANDOCK_VERSION}/plugins
)

