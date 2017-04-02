project(candock)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

set( CANDOCK_MAJOR_VERSION 0 )
set( CANDOCK_MINOR_VERSION 3 )
set( CANDOCK_TWEAK_VERSION 3 )

set(CANDOCK_VERSION "v${CANDOCK_MAJOR_VERSION}.${CANDOCK_MINOR_VERSION}.${CANDOCK_TWEAK_VERSION}"
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

#Use included version if not installed on compile machine
if ( NOT RAPIDXML_ROOT )
    set( RAPIDXML_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/foreign/ )
endif( NOT RAPIDXML_ROOT )

find_package(Boost    REQUIRED  COMPONENTS date_time filesystem system regex program_options)
find_package(GSL      REQUIRED)
find_package(RapidXML REQUIRED)
find_package(OpenMM   REQUIRED)
find_package(Git      REQUIRED)

set(CANDOCK_CUDA_SUPPORT OFF CACHE BOOL "Build CUDA support")

if (CANDOCK_CUDA_SUPPORT)

    find_package(CUDA REQUIRED)

    if(CUDA_FOUND)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=compute_61 -code=sm_61" )
        add_definitions(-DCOMPILE_CUDA)
    endif(CUDA_FOUND)

endif(CANDOCK_CUDA_SUPPORT)

set(JSONCPP_INCLUDE_DIR "foreign/jsoncpp/include" "foreign/jsoncpp/src" )

if (WIN32)

    set(Boost_USE_STATIC_LIBS        OFF)
    set(Boost_USE_MULTITHREADED      ON)    
    set(Boost_USE_STATIC_RUNTIME     OFF)
    set(BOOST_ALL_DYN_LINK           ON)   # force dynamic linking for all libraries

    add_definitions( /D_USE_MATH_DEFINES /DBOOST_ALL_NO_LIB /DBOOST_ALL_DYN_LINK /DWIN32_LEAN_AND_MEAN)

    if (NOT CANDOCK_BUILD_STATIC_EXECUTABLE)
	    add_definitions( /DCANDOCK_SHARED_LIBRARY /DKB_USE_STATIC_LIBRARIES)
    endif(NOT CANDOCK_BUILD_STATIC_EXECUTABLE)

else(WIN32)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -fPIC")

endif(WIN32)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIR})
include_directories(${RAPIDXML_DIR})
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

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(scripts)

set(CANDOCK_BUILD_GUI OFF CACHE BOOL "Build the CANDOCK gui")

if (CANDOCK_BUILD_GUI)
    add_subdirectory(gui)
endif (CANDOCK_BUILD_GUI)

add_subdirectory(doc)
