project(candock)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

set( CANDOCK_MAJOR_VERSION 0 )
set( CANDOCK_MINOR_VERSION 3 )
set( CANDOCK_TWEAK_VERSION 0 )

set(CANDOCK_VERSION "v${CANDOCK_MAJOR_VERSION}.${CANDOCK_MINOR_VERSION}.${CANDOCK_BUILD_VERSION}"
     CACHE STRING "CANDOCK version." FORCE
)

mark_as_advanced(CANDOCK_VERSION)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

set(CANDOCK_BUILD_STATIC_EXECUTABLE OFF CACHE BOOL "Build a version of candock without dynamic linking")

if(${CANDOCK_BUILD_STATIC_EXECUTABLE})
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif(WIN32)
endif(${CANDOCK_BUILD_STATIC_EXECUTABLE})

find_package(Boost   REQUIRED  COMPONENTS date_time filesystem system regex program_options)
find_package(GSL     REQUIRED)
find_package(OpenMM  REQUIRED)
find_package(Git     REQUIRED)
find_package(CUDA)

if(CUDA_FOUND)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=compute_61 -code=sm_61" )
    add_definitions(-DCOMPILE_CUDA)
endif(CUDA_FOUND)


#Change to c++11 due to gnu specific extensions that break cuda
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -fPIC")
set(JSONCPP_INCLUDE_DIR "foreign/jsoncpp/include" "foreign/jsoncpp/src" )

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIR})
include_directories(${OPENMM_INCLUDE_DIR})
include_directories(lib)

execute_process(
    COMMAND ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
    OUTPUT_VARIABLE CANDOCK_GIT_REFERENCE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE CANDOCK_GIT_MYCBRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    ${PROJECT_SOURCE_DIR}/src/version.hpp.in
    ${PROJECT_BINARY_DIR}/version.hpp
)

include_directories( ${PROJECT_BINARY_DIR} )

include_directories( foreign )

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(scripts)
add_subdirectory(doc)
