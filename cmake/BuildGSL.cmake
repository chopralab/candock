include(ExternalProject)

set(USING_INTERNAL_GSL on CACHE BOOL "Using Internal GSL")
mark_as_advanced(USING_INTERNAL_GSL)

ExternalProject_Add( GSL
    URL https://github.com/ampl/gsl/archive/644e768630841bd085cb7121085a688c4ff424d0.tar.gz
    URL_HASH SHA1=de333ea777c1880edb9a59d68c8ed3e8cc91810c
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsl
    CMAKE_CACHE_ARGS    -DGSL_DISABLE_TESTS:Bool=on
                        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/gsl_build
                        -DGSL_DISABLE_WARNINGS:Bool=on
                        -DBUILD_SHARED_LIBS:Bool=on
)

add_library(gsl SHARED IMPORTED)
add_library(gslcblas SHARED IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gsl${CMAKE_SHARED_LIBRARY_SUFFIX})
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/libgslcblas.so)
add_dependencies(gsl GSL)
add_dependencies(gslcblas GSL)

set(GSL_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/include CACHE PATH "Included GSL Dir")
link_directories(${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib)
set(GSL_LIBRARIES gsl gslcblas)

mark_as_advanced(GSL_INCLUDE_DIRS)
mark_as_advanced(GSL_LIBRARIES)

install(FILES
 ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gsl${CMAKE_SHARED_LIBRARY_SUFFIX}
 ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gslcblas${CMAKE_SHARED_LIBRARY_SUFFIX}
DESTINATION
 ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/lib
)
