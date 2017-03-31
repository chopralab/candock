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

find_package(Boost   REQUIRED  COMPONENTS date_time filesystem system regex program_options)
find_package(GSL     REQUIRED)
find_package(OpenMM  REQUIRED)
find_package(Git     REQUIRED)

set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
set(BOOST_ALL_DYN_LINK           ON)   # force dynamic linking for all libraries

#Change to c++11 due to gnu specific extensions that break cuda
set(JSONCPP_INCLUDE_DIR "foreign/jsoncpp/include" "foreign/jsoncpp/src" )
add_definitions( /D__USE_MATH_DEFINES /DBOOST_ALL_NO_LIB /DBOOST_ALL_DYN_LINK)

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
