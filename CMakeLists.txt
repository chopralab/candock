project(candock)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

find_package(Boost   REQUIRED  COMPONENTS date_time filesystem system regex)
find_package(GSL     REQUIRED)
find_package(OpenMM REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(JSONCPP_INCLUDE_DIR "foreign/jsoncpp-src-0.6.0-rc2/include" "foreign/jsoncpp-src-0.6.0-rc2/src" )
set(RAPIDXML_INCLUDE_DIR "foreign/" )

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIR})
include_directories(${RAPIDXML_INCLUDE_DIR})

include_directories(lib)

add_definitions( -DNDEBUG )

include_directories(${OPENMM_INCLUDE_DIR})
include_directories(lib)

#~ link_directories( ${CMAKE_CURRENT_SOURCE_DIR}/openmm/lib )

add_subdirectory(lib)
add_subdirectory(src)
