#---------------------------------------------------
# OpenMM Example Plugin
#----------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

set(KBFORCE_CUDA_SUPPORT OFF CACHE BOOL "Build CUDA support")

IF(KBFORCE_CUDA_SUPPORT)
    FIND_PACKAGE(CUDA)
ENDIF(KBFORCE_CUDA_SUPPORT)

IF(CUDA_FOUND AND KBFORCE_CUDA_SUPPORT)
    SET(EXAMPLE_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND AND KBFORCE_CUDA_SUPPORT)
    SET(EXAMPLE_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND AND KBFORCE_CUDA_SUPPORT)

INCLUDE_DIRECTORIES("${OPENMM_INCLUDE_DIR}")
link_directories("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# set flags for linking on mac
IF(APPLE)
    SET (CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-msse2 -stdlib=libc++")
ENDIF(APPLE)

IF(WIN32)
    ADD_DEFINITIONS( /DKB_BUILDING_SHARED_LIBRARY )
ENDIF(WIN32)

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX ${OPENMM_DIR} CACHE PATH "Where to install the plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(EXAMPLE_PLUGIN_SOURCE_SUBDIRS openmmapi serialization)

# Set the library name
SET(EXAMPLE_LIBRARY_NAME KBPlugin)
SET(SHARED_EXAMPLE_TARGET ${EXAMPLE_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${EXAMPLE_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Create the library.

ADD_LIBRARY(${SHARED_EXAMPLE_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_EXAMPLE_TARGET}
    PROPERTIES COMPILE_FLAGS "-DEXAMPLE_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_EXAMPLE_TARGET} ${OPENMM_LIBRARY})

# install headers
#FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
#INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
#FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
#INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

# Enable testing

#ENABLE_TESTING()
#ADD_SUBDIRECTORY(serialization/tests)

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

IF(EXAMPLE_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(EXAMPLE_BUILD_CUDA_LIB)

