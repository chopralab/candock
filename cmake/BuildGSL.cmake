include(ExternalProject)
ExternalProject_Add( GSL
    URL ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
    URL_HASH SHA1=5db17d71306139b811a218d8a8cf98e12e1f58ae
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsl
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gsl/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/gsl_build
    BUILD_COMMAND ${MAKE}
)

set(GSL_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/include CACHE PATH "Included GSL Dir")
link_directories(${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib)
set(GSL_LIBRARIES gsl gslcblas)

install(FILES
 ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/libgsl.so.23
 ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/libgsl.so.23.0.0
 ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/libgslcblas.so.0
 ${CMAKE_CURRENT_BINARY_DIR}/gsl_build/lib/libgslcblas.so.0.0.0
DESTINATION
 ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/lib
)