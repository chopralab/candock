# Update this value if you need to update the data file set

set(TESTS_DATA_GIT "8fe2b5a5c4c399cdebdf9e6950c2ead9a17b586f")

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/external/openmm.tar.gz")
    message(STATUS "Downloading test data files")
    file(DOWNLOAD
        "https://github.com/pandegroup/openmm/archive/${TESTS_DATA_GIT}.tar.gz"
        "${CMAKE_CURRENT_BINARY_DIR}/openmm.tar.gz"
        SHOW_PROGRESS
        EXPECTED_HASH SHA1=4c7c70e6113d31e11a9efe574e6626619e0de34b
    )

    message(STATUS "Unpacking test data files")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E remove_directory openmm
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xf openmm.tar.gz
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )


    execute_process(
        COMMAND ${CMAKE_COMMAND} -E rename openmm-${TESTS_DATA_GIT} openmm
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E touch openmm/${TESTS_DATA_GIT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endif()

set(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS OFF CACHE BOOL "Unwanted OpenMM Feature")
set(OPENMM_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Unwanted OpenMM Feature")
set(OPENMM_BUILD_EXAMPLES OFF CACHE BOOL "Unwanted OpenMM Feature")
set(OPENMM_BUILD_REFERENCE_TESTS OFF CACHE BOOL "Unwanted OpenMM Feature")

set(OPENMM_BUILD_SHARED_LIB ON CACHE BOOL "Do not build OpenMM as a shared library")
set(OPENMM_BUILD_STATIC_LIB OFF CACHE BOOL "DO build it statically")

SET(OPENMM_BUILD_AMOEBA_PLUGIN OFF CACHE BOOL "Build Amoeba plugin")
SET(OPENMM_BUILD_RPMD_PLUGIN ON CACHE BOOL "Build RPMD plugin")
SET(OPENMM_BUILD_DRUDE_PLUGIN ON CACHE BOOL "Build Drude plugin")
SET(OPENMM_BUILD_PME_PLUGIN OFF CACHE BOOL "Build CPU PME plugin")
SET(OPENMM_BUILD_CUDA_COMPILER_PLUGIN OFF CACHE BOOL "Build CUDA runtime compiler plugin")

add_definitions(-DIEEE_8087)

add_subdirectory(
    ${CMAKE_CURRENT_BINARY_DIR}/openmm
    EXCLUDE_FROM_ALL
)

set(OPENMM_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmm/openmmapi/include)
list(APPEND OPENMM_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmm/olla/include)
list(APPEND OPENMM_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmm/platforms/reference/include)
list(APPEND OPENMM_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/foreign/kbforce/openmmapi/include)
set(OPENMM_LIBRARY OpenMM CACHE STRING "Internal OpenMM Library" )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/foreign/kbforce)

set(OPENMM_KBPLUGIN KBPlugin CACHE STRING  "KBForce Internal Library")

mark_as_advanced(OPENMM_LIBRARY)

install(FILES
 ${CMAKE_CURRENT_BINARY_DIR}/libOpenMM.so
DESTINATION
 ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/lib
)
