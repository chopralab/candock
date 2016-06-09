project(candock)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

find_package(Boost    REQUIRED  COMPONENTS date_time filesystem regex)
find_package(GSL      REQUIRED)
find_package(Jsoncpp  REQUIRED)
find_package(RapidXML REQUIRED)
find_package(OpenMM   REQUIRED)
find_package(TCLAP    REQUIRED)

set (CMAKE_CXX_STANDARD 11)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIR})
include_directories(${RAPIDXML_DIR})
include_directories(${OpenMM_INCLUDE_DIR})
include_directories(${TCLAP_INCLUDE_DIR})

include_directories(lib)

add_definitions( -DNDEBUG )

add_subdirectory(lib)
add_subdirectory(src)
