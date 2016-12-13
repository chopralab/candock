project(candock)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set( CANDOCK_MAJOR_VERSION 0 )
set( CANDOCK_MINOR_VERSION 2 )
set( CANDOCK_TWEAK_VERSION 1 )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

IF(WIN32)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
ELSE(WIN32)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
ENDIF(WIN32)

find_package(Boost   REQUIRED  COMPONENTS date_time filesystem system regex program_options)
find_package(GSL     REQUIRED)
find_package(OpenMM  REQUIRED)
find_package(Git     REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
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

configure_file( ${PROJECT_SOURCE_DIR}/src/version.hpp.in
                ${PROJECT_BINARY_DIR}/version.h
)
            
include_directories( ${PROJECT_BINARY_DIR} )

include_directories( foreign )

add_subdirectory(lib)
add_subdirectory(src)
